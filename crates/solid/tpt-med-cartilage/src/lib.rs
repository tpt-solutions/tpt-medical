//! Biphasic / poroelastic cartilage models.
//!
//! Linear biphasic theory (Mow, Kuei & Lai 1980) in the confined-compression
//! configuration: a solid matrix (aggregate modulus `H_A`, Poisson ratio
//! ≈ 0) saturated by interstitial fluid, with Darcy permeability `k`.
//!
//! Creep under a step load `σ0` has the closed-form surface displacement
//!
//! ```text
//! u(t)/h = 1 − (1 − δ(t))·(1 − δ(0))...
//! u(t)/h = (σ0/H_A) · [1 − (1+2h²/(π²kt))⁻¹·...]
//! ```
//!
//! implemented here via the classical series solution; the short-time
//! (fluid-supported, t→0) and long-time (solid-supported, t→∞) limits are
//! verified against the analytic limits `u/h → σ0/H_A` and `u(0)=0`.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

/// A strain-dependent permeability law: `k` as a function of a
/// deformation measure. The classical cartilage physics — permeability
/// falling steeply as the matrix compacts — is what makes biphasic
/// load responses nonlinear in time even when the solid matrix is linear.
///
/// The crate ships the constant-`k` default (the linear biphasic
/// baseline); strain-dependent forms (exponential in compaction, or in
/// `J²`, per the literature's several parameterisations) are
/// caller-supplied closures implementing this trait, with their own
/// citations — the same pattern as `tpt-med-bone`'s `ModulusLaw` and
/// `tpt-med-dicom`'s calibration hooks: the mechanism ships, the
/// coefficients come from the caller's cited source.
pub trait PermeabilityLaw {
    /// Permeability `k` (mm⁴/(N·s)) at a volume ratio `J` (dimensionless;
    /// `J < 1` = compaction).
    fn permeability(&self, volume_ratio: f64) -> f64;
}

/// The linear-biphasic baseline: constant [`BiphasicMaterial::permeability`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConstantPermeability {
    /// The constant `k` (mm⁴/(N·s)).
    pub k: f64,
}

impl PermeabilityLaw for ConstantPermeability {
    fn permeability(&self, _volume_ratio: f64) -> f64 {
        self.k
    }
}

impl<F: Fn(f64) -> f64> PermeabilityLaw for F {
    fn permeability(&self, volume_ratio: f64) -> f64 {
        self(volume_ratio)
    }
}

/// A biphasic material with a caller-supplied strain-dependent
/// permeability, for creep/relaxation sweeps under matrix compaction.
///
/// The evaluation points this crate owns (the creep series and its
/// `gel_time`) are closed-form solutions of the CONSTANT-`k` problem; a
/// strain-dependent `k` makes the diffusion coefficient time- and
/// space-dependent and those series no longer apply. What this type
/// therefore provides is the **effective-permeability evaluation**: the
/// `k` at a given compaction level, for callers driving their own
/// stepping scheme zone by zone.
#[derive(Debug, Clone)]
pub struct StrainDependentPermeability<L: PermeabilityLaw> {
    /// The constant-`k` baseline (also the reference the law should be
    /// calibrated against: `law.permeability(1.0)` should be the
    /// unconstrained-compaction value).
    pub material: BiphasicMaterial,
    /// The permeability law.
    pub law: L,
}

impl<L: PermeabilityLaw> StrainDependentPermeability<L> {
    /// Permeability at a volume ratio, validated positive and finite.
    pub fn permeability_at(&self, volume_ratio: f64) -> Result<f64, String> {
        let k = self.law.permeability(volume_ratio);
        if !k.is_finite() || k <= 0.0 {
            return Err(format!(
                "permeability law returned {k} at J = {volume_ratio}"
            ));
        }
        Ok(k)
    }

    /// Permeability at the equilibrium compaction of a step stress `σ0`
    /// under the constant-`k` aggregate-modulus relation — the screening
    /// point at which a creep sweep's late-time `k` is evaluated.
    pub fn permeability_at_equilibrium(&self, sigma0: f64) -> Result<f64, String> {
        let j = 1.0 - self.material.equilibrium_strain(sigma0);
        self.permeability_at(j)
    }
}

/// Biphasic cartilage material.
#[derive(Debug, Clone, Copy)]
pub struct BiphasicMaterial {
    /// Aggregate modulus `H_A` (MPa).
    pub aggregate_modulus: f64,
    /// Darcy permeability `k` (mm⁴/(N·s)).
    pub permeability: f64,
    /// Solid-matrix Poisson ratio (≈ 0 for cartilage).
    pub poissons_ratio: f64,
    /// Sample thickness in the loading direction (mm).
    pub thickness: f64,
}

impl Default for BiphasicMaterial {
    fn default() -> Self {
        // Screening values for adult articular cartilage.
        Self {
            aggregate_modulus: 0.7,
            permeability: 0.002,
            poissons_ratio: 0.0,
            thickness: 2.0,
        }
    }
}

impl BiphasicMaterial {
    /// Confined-compression creep series: surface displacement fraction
    /// `u(t)/h` under a step stress `σ0` (MPa), classical solution
    ///
    /// ```text
    /// u(t)/h = (σ0/H_A) [ 1 − (8/π²) Σ_{n odd} (1/n²) exp(−n²π² k H_A t / (4 h²)) ]
    /// ```
    ///
    /// evaluated with `terms` series terms (16 is converged to machine
    /// precision for practical t).
    pub fn creep_displacement_fraction(&self, sigma0: f64, time: f64, terms: usize) -> f64 {
        let h = self.thickness;
        let diffusivity = self.permeability * self.aggregate_modulus / (4.0 * h * h);
        let mut series = 0.0;
        for n in 0..terms {
            let m = (2 * n + 1) as f64;
            series += (1.0 / (m * m))
                * (-m * m * core::f64::consts::PI * core::f64::consts::PI * diffusivity * time)
                    .exp();
        }
        (sigma0 / self.aggregate_modulus)
            * (1.0 - 8.0 / (core::f64::consts::PI * core::f64::consts::PI) * series)
    }

    /// Half-equilibrium ("gel") time constant
    /// `t½ = 4 ln2 · h² / (π² k H_A)` (s): the confined-compression creep
    /// reaches half of its equilibrium displacement at this time
    /// (single-mode approximation of the series solution).
    pub fn gel_time(&self) -> f64 {
        4.0 * (2.0f64).ln() * self.thickness * self.thickness
            / (core::f64::consts::PI
                * core::f64::consts::PI
                * self.permeability
                * self.aggregate_modulus)
    }

    /// Instantaneous (fluid-supported) modulus at `t → 0+`: the response is
    /// rigid; displacement starts at zero.
    pub fn initial_displacement_fraction(&self, sigma0: f64) -> f64 {
        self.creep_displacement_fraction(sigma0, 0.0, 32)
    }

    /// Equilibrium strain `σ0/H_A` at `t → ∞`.
    pub fn equilibrium_strain(&self, sigma0: f64) -> f64 {
        sigma0 / self.aggregate_modulus
    }

    /// Unconfined-compression **equilibrium** modulus (MPa): the stress
    /// at `t → ∞` per unit strain when the free-draining side has
    /// depressurised and the solid matrix carries everything alone —
    /// `E_s = H_A·(1+ν_s)(1−2ν_s)/(1−ν_s)`, the solid matrix's Young's
    /// modulus. For the cartilage default `ν_s = 0` this equals `H_A`
    /// exactly, so the confined and unconfined equilibria coincide.
    /// `None` for `ν_s ≥ ½` (no valid linear modulus).
    ///
    /// The unconfined **transient** between the rigid instantaneous
    /// response (fluid-supported, `u(0) = 0`) and this equilibrium is the
    /// classical Bessel-series solution of the coupled radial/axial
    /// problem and is deliberately not implemented here — its
    /// coefficients are cited literature, not derivations to be
    /// reproduced from memory. The two limits are exact and testable;
    /// the transient needs the series.
    pub fn unconfined_equilibrium_modulus(&self) -> Option<f64> {
        let nu = self.poissons_ratio;
        if nu < 0.0 || nu >= 0.5 {
            return None;
        }
        Some(self.aggregate_modulus * (1.0 + nu) * (1.0 - 2.0 * nu) / (1.0 - nu))
    }

    /// Unconfined-compression **instantaneous** modulus (MPa): the
    /// response at `t → 0⁺`, before any interstitial flow — "the biphasic
    /// continuum deforms without change in volume and behaves like an
    /// incompressible elastic solid of the same shear modulus"
    /// (Armstrong, Lai & Mow 1984, *J Biomech Eng* 106:165, abstract). An
    /// incompressible isotropic solid with shear modulus `G` has
    /// `E = 3G`, so the bookends of the unconfined transient are
    /// `E(0⁺) = 3G` here and [`Self::unconfined_equilibrium_modulus`] =
    /// `2G(1+ν_s)` at `t → ∞`: at the cartilage default `ν_s = 0` the
    /// wall relaxes from `1.5·H_A` down to `H_A`, and the relaxation
    /// ratio `3/(2(1+ν_s))` diverges as `ν_s → ½` (the unconfined
    /// response becomes confined-like), which is the classical
    /// sensitivity the transient series exists to resolve in time. The
    /// transient itself stays deferred — its Bessel-series coefficients
    /// are cited literature (paywalled), not memory-reproducible
    /// derivations. `None` for `ν_s ≥ ½`.
    pub fn unconfined_instantaneous_modulus(&self) -> Option<f64> {
        self.solid_shear_modulus().map(|g| 3.0 * g)
    }

    /// Solid-matrix shear modulus `G = H_A·(1−2ν_s)/(2(1−ν_s))` (MPa).
    ///
    /// First-order biphasic **shear carries no interstitial fluid
    /// pressurisation**: shear produces no volumetric strain, so Darcy
    /// flow has nothing to drive and the response is the solid matrix at
    /// *every* time — unlike compression, there is no transient and no
    /// boundary-condition choice to make. The shear half of the
    /// unconfined-shear boundary-condition item is therefore a closed
    /// form, not a solve. `None` for `ν_s ≥ ½`, where the linear theory
    /// has no valid shear modulus.
    pub fn solid_shear_modulus(&self) -> Option<f64> {
        let nu = self.poissons_ratio;
        if nu < 0.0 || nu >= 0.5 {
            return None;
        }
        Some(self.aggregate_modulus * (1.0 - 2.0 * nu) / (2.0 * (1.0 - nu)))
    }

    /// Interstitial fluid pressure fraction `p(t)/σ0` at the impermeable
    /// subchondral boundary:
    /// `p/σ0 = (8/π) Σ_{n odd} (1/n) exp(−n² t/t_c)·sin(nπ·x/h)` evaluated
    /// at the sealed surface `x = 0` (sin → 0) — here evaluated at mid-plane
    /// `x = h/2`, the standard screening location.
    pub fn fluid_pressure_fraction(&self, time: f64, terms: usize) -> f64 {
        let h = self.thickness;
        let diffusivity = self.permeability * self.aggregate_modulus / (4.0 * h * h);
        let mut series = 0.0;
        for n in 0..terms {
            let m = (2 * n + 1) as f64;
            series += (1.0 / m)
                * (-m * m * core::f64::consts::PI * core::f64::consts::PI * diffusivity * time)
                    .exp()
                * (m * core::f64::consts::PI * 0.5).sin();
        }
        8.0 / core::f64::consts::PI * series
    }
}

/// Squeeze-film lubrication of the contact interface: a Newtonian film of
/// viscosity `μ` squeezed between parallel circular bearing surfaces of
/// radius `R` — the screening geometry for synovial-joint lubrication,
/// where the interstitial fluid phase is the lubricant and the film
/// carries the load while it is being squeezed out.
///
/// The constitutive relation is **Stefan's equation** (1874), the
/// parallel-plate limit of Reynolds' equation: with the pressure solving
/// `∇²p = 12μ|ḣ|/h³` over the disk (`p(R) = 0`), the load the film
/// carries at thickness `h` under approach rate `−ḣ` is
///
/// ```text
/// W = (3π μ R⁴ / 2 h³) · (−ḣ)
/// ```
///
/// and a constant load `W` therefore squeezes the film as
/// `h(t) = [h₀⁻² + 4Wt/(3πμR⁴)]^(−1/2)` — fast at first, ever slower,
/// never touching (the classical squeeze-film cushion). This is the
/// lubrication term the contact item asked for: [`Self::load_capacity`]
/// is the constitutive term a contact solver evaluates (approach rate in,
/// carried load out), and [`Self::film_thickness`] / [`Self::time_to_squeeze`]
/// are the step-load creep forms. Viscosity is caller-supplied (synovial
/// fluid is strongly shear-rate dependent; a screening constant in the
/// 0.005–0.5 Pa·s band converts to 5e-9–5e-7 MPa·s in this crate's
/// N-mm-s units) — the mechanism ships, the coefficient comes from the
/// caller's source.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SqueezeFilm {
    /// Lubricant dynamic viscosity `μ` (MPa·s = N·s/mm²).
    pub viscosity: f64,
    /// Bearing-surface contact radius `R` (mm).
    pub radius: f64,
}

impl SqueezeFilm {
    /// Validates; `Err` for a non-finite or non-positive viscosity or
    /// radius.
    pub fn new(viscosity: f64, radius: f64) -> Result<Self, String> {
        if !viscosity.is_finite() || viscosity <= 0.0 {
            return Err(format!(
                "viscosity must be positive and finite, got {viscosity}"
            ));
        }
        if !radius.is_finite() || radius <= 0.0 {
            return Err(format!("radius must be positive and finite, got {radius}"));
        }
        Ok(Self { viscosity, radius })
    }

    /// Film thickness `h(t)` (mm) after `time` (s) under a constant load
    /// `load` (N), from an initial film thickness (mm). Stefan's equation;
    /// `h(0) = initial`, monotonically decreasing, never zero.
    pub fn film_thickness(&self, initial: f64, load: f64, time: f64) -> f64 {
        let inv_sq = initial * initial;
        (1.0 / inv_sq
            + 4.0 * load * time
                / (3.0 * core::f64::consts::PI * self.viscosity * self.radius.powi(4)))
        .sqrt()
        .recip()
    }

    /// The load (N) the film carries at `thickness` (mm) under approach
    /// rate `approach_rate` = `−dh/dt` (mm/s) — the contact-solver-facing
    /// constitutive term.
    pub fn load_capacity(&self, thickness: f64, approach_rate: f64) -> f64 {
        3.0 * core::f64::consts::PI * self.viscosity * self.radius.powi(4)
            / (2.0 * thickness.powi(3))
            * approach_rate
    }

    /// Time (s) to squeeze the film under constant `load` from one
    /// thickness to a thinner one (mm) — the exact integral of Stefan's
    /// equation. `Err` when `to ≥ from` (no thinning) or either is
    /// non-positive.
    pub fn time_to_squeeze(&self, load: f64, from: f64, to: f64) -> Result<f64, String> {
        if from <= 0.0 || to <= 0.0 {
            return Err(format!("thicknesses must be positive, got {from} -> {to}"));
        }
        if to >= from {
            return Err(format!("the film thins: {from} -> {to} is no squeeze"));
        }
        if !load.is_finite() || load <= 0.0 {
            return Err(format!("load must be positive and finite, got {load}"));
        }
        Ok(
            3.0 * core::f64::consts::PI * self.viscosity * self.radius.powi(4) / (4.0 * load)
                * (1.0 / (to * to) - 1.0 / (from * from)),
        )
    }
}

/// A 1-D confined-compression creep stepper for an arbitrary
/// [`PermeabilityLaw`] — the time-stepping half of the nonlinear-biphasic
/// item. The closed-form series in [`BiphasicMaterial`]
/// (`creep_displacement_fraction`, `fluid_pressure_fraction`) are
/// solutions of the CONSTANT-`k` problem; a strain-dependent `k` makes
/// the diffusion coefficient space- and time-dependent and those series
/// no longer apply. This stepper integrates the same problem numerically:
///
/// ```text
/// ṗ = k(J) H_A ∂²p/∂x²,   e = (σ0 − p)/H_A,   J = 1 − e
/// ```
///
/// on `cells` nodes over the layer (drained surface at `x = 0`, sealed
/// platen at `x = h` — the same boundary convention the series uses), so
/// `u(t)/h = (σ0 − p̄)/H_A` throughout. When the law is
/// [`ConstantPermeability`], the stepper reproduces the closed-form
/// series to discretisation error — which is the verification that makes
/// the strain-dependent results trustworthy. Integration is explicit
/// Euler with an adaptive step held to a fraction of the diffusion CFL
/// bound computed from the current `k` field (no monotonicity assumption
/// on the law).
#[derive(Debug, Clone)]
pub struct ConfinedCreepStepper {
    /// The (linear-elastic solid matrix) biphasic material.
    pub material: BiphasicMaterial,
    /// Spatial resolution: number of cells over the layer thickness.
    pub cells: usize,
}

impl ConfinedCreepStepper {
    /// Validates; `Err` for fewer than 8 cells (the boundary layer the
    /// series resolves is unresolvable coarser than that).
    pub fn new(material: BiphasicMaterial, cells: usize) -> Result<Self, String> {
        if cells < 8 {
            return Err(format!("at least 8 cells are needed, got {cells}"));
        }
        Ok(Self { material, cells })
    }

    fn cell_permeabilities(&self, p: &[f64], sigma0: f64, law: &impl PermeabilityLaw) -> Vec<f64> {
        let ha = self.material.aggregate_modulus;
        p.iter()
            .map(|&p_i| {
                let j = 1.0 - (sigma0 - p_i) / ha;
                law.permeability(j)
            })
            .collect()
    }

    /// One explicit-Euler step of `p` (length `cells + 1`; `p[0]` is the
    /// drained surface and stays zero). `dt` must satisfy the
    /// [`Self::stable_time_step`] bound.
    pub fn step(&self, p: &mut [f64], sigma0: f64, dt: f64, law: &impl PermeabilityLaw) {
        let n = self.cells;
        debug_assert_eq!(p.len(), n + 1, "one pressure value per node");
        let ha = self.material.aggregate_modulus;
        let k = self.cell_permeabilities(p, sigma0, law);
        let dx2 = (self.material.thickness / n as f64).powi(2);
        let mut dp = vec![0.0; n + 1];
        for i in 1..n {
            dp[i] = ha * k[i] * (p[i + 1] - 2.0 * p[i] + p[i - 1]) / dx2;
        }
        // Sealed platen: mirrored node, so the second difference collapses
        // to (p[n-1] - p[n]).
        dp[n] = ha * k[n] * (p[n - 1] - p[n]) / dx2;
        for i in 1..=n {
            p[i] += dt * dp[i];
        }
    }

    /// The diffusion CFL bound for the current state: the explicit step
    /// must stay below `dx² / (2 · max(k) · H_A)`; the driver uses 0.4 of
    /// it.
    pub fn stable_time_step(&self, p: &[f64], sigma0: f64, law: &impl PermeabilityLaw) -> f64 {
        let k = self.cell_permeabilities(p, sigma0, law);
        let k_max = k.iter().cloned().fold(0.0f64, f64::max);
        let dx = self.material.thickness / self.cells as f64;
        dx * dx / (2.0 * k_max * self.material.aggregate_modulus)
    }

    /// The pressure field after integration to `time` — for callers who
    /// want the profile, not just the scalar creep fraction. The field is
    /// `[p(0) = 0, p_1, ..., p(cells)]` over the layer.
    pub fn pressure_profile_with_law(
        &self,
        sigma0: f64,
        time: f64,
        law: &impl PermeabilityLaw,
    ) -> Vec<f64> {
        let mut p = vec![sigma0; self.cells + 1];
        p[0] = 0.0;
        let mut t = 0.0f64;
        while t < time {
            let dt = (0.8 * self.stable_time_step(&p, sigma0, law)).min(time - t);
            self.step(&mut p, sigma0, dt, law);
            p[0] = 0.0;
            t += dt;
        }
        p
    }

    /// Mean strain `u(t)/h` of a pressure field — the creep displacement
    /// fraction (zero at `t = 0`, `σ0/H_A` at equilibrium).
    pub fn creep_fraction_of(&self, p: &[f64], sigma0: f64) -> f64 {
        let ha = self.material.aggregate_modulus;
        let p_bar = p.iter().sum::<f64>() / p.len() as f64;
        (sigma0 - p_bar) / ha
    }

    /// Integrates confined-compression creep under `sigma0` (MPa) to
    /// `time` (s) with the given permeability law, returning `u(t)/h`.
    /// The step is adaptive (a fixed fraction of the state-dependent CFL
    /// bound), so the call is deterministic for a given `cells`.
    pub fn creep_fraction_with_law(
        &self,
        sigma0: f64,
        time: f64,
        law: &impl PermeabilityLaw,
    ) -> f64 {
        let p = self.pressure_profile_with_law(sigma0, time, law);
        self.creep_fraction_of(&p, sigma0)
    }
}

/// A 1-D confined-compression stepper for a **nonlinear solid matrix**:
/// the drained effective stress `σ_eff(e)` is a caller-supplied
/// increasing function of the apparent compressive strain
/// `e = 1 − J`, replacing the linear `σ_eff = H_A·e` the
/// [`ConfinedCreepStepper`] embeds.
///
/// # Formulation
///
/// The quasi-static confined problem has one independent field once the
/// equilibrium `σ_eff(e) + p = σ₀` is substituted (effective solid
/// stress plus excess pore pressure equals the applied `σ₀`; undrained
/// `t = 0⁺` has `e = 0`, `p = σ₀` — the fluid carries everything).
/// Darcy continuity then gives a heat-like equation for the excess
/// pressure with the **state-dependent diffusivity**
///
/// ```text
/// ∂p/∂t = k(e) · σ_eff′(e) · p_xx,    e = σ_eff⁻¹(σ₀ − p)
/// ```
///
/// the poroelastic coefficient that generalizes the linear `k·H_A`. The
/// drained surface holds the constant Dirichlet pressure `p = 0` (which
/// pins `e(0) = σ_eff⁻¹(σ₀)` = the equilibrium strain); the platen end
/// is sealed (`p_x = 0`). The creep fraction is `u/h = mean(e)` — zero
/// undrained, `σ_eff⁻¹(σ₀)` at equilibrium.
///
/// `σ_eff` is caller-supplied and caller-cited (the mechanism ships,
/// the coefficients come from the caller's source, as everywhere in
/// this workspace); it must be increasing with `σ_eff(0) = 0`. The
/// update is explicit with an adaptive step bounded by the
/// state-dependent diffusion limit `dx²/(2·max k·σ_eff′)`.
#[derive(Debug, Clone)]
pub struct NonlinearConfinedStepper {
    /// Tissue thickness in the loading direction (mm).
    pub thickness: f64,
    /// Spatial resolution: cells over the thickness.
    pub cells: usize,
}

/// Inverts an increasing `σ_eff` on the compressive-strain range by a
/// bracketed bisection whose bracket expands from zero.
fn invert_sigma_eff(sigma_eff: &impl Fn(f64) -> f64, target: f64) -> Result<f64, String> {
    if target < 0.0 {
        return Err(format!("negative driving stress {target}"));
    }
    let mut lo = 0.0f64;
    let mut hi = 0.05f64;
    while sigma_eff(hi) < target {
        lo = hi;
        hi *= 2.0;
        if hi > 4.0 {
            return Err(format!(
                "sigma_eff does not reach {target} by e = 4 (300 % compressive strain)"
            ));
        }
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if sigma_eff(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

impl NonlinearConfinedStepper {
    /// Validates; `Err` for fewer than 8 cells.
    pub fn new(thickness: f64, cells: usize) -> Result<Self, String> {
        if cells < 8 {
            return Err(format!("at least 8 cells are needed, got {cells}"));
        }
        Ok(Self { thickness, cells })
    }

    /// `σ_eff′` by central differences (the caller's law is a black box;
    /// the step is scaled by `max(1, |e|)` like the tissue crate's
    /// tangent differencing).
    fn sigma_eff_prime(sigma_eff: &impl Fn(f64) -> f64, e: f64, h: f64) -> f64 {
        let step = h * e.abs().max(1.0);
        (sigma_eff(e + step) - sigma_eff(e - step)) / (2.0 * step)
    }

    /// The strain field `e(x)` after marching to `time` under `σ0` with
    /// the given permeability and drained-stress laws. The drained
    /// surface holds `e(0) = σ_eff⁻¹(σ0)` from `t > 0`.
    pub fn strain_profile(
        &self,
        sigma0: f64,
        time: f64,
        permeability: &impl PermeabilityLaw,
        sigma_eff: &impl Fn(f64) -> f64,
    ) -> Result<Vec<f64>, String> {
        let n = self.cells;
        let dx = self.thickness / n as f64;
        // The drained surface's equilibrium strain; asserted reachable so
        // an unreachable sigma0 fails before the march, not during it.
        let _e_eq = invert_sigma_eff(sigma_eff, sigma0)?;

        // Undrained start: zero strain, excess pressure sigma0 everywhere,
        // drained surface pinned at p = 0.
        let mut p = vec![sigma0; n + 1];
        p[0] = 0.0;
        let mut t = 0.0f64;
        while t < time {
            // Per-cell strain and state-dependent diffusivity
            // D = k(e)·sigma_eff'(e).
            let mut d_node = vec![0.0f64; n + 1];
            for i in 1..=n {
                let e_i = invert_sigma_eff(sigma_eff, sigma0 - p[i])?;
                let prime = Self::sigma_eff_prime(sigma_eff, e_i, 1.0e-6);
                if !prime.is_finite() || prime <= 0.0 {
                    return Err(format!(
                        "sigma_eff must be increasing: sigma_eff'({e_i:.4}) = {prime}"
                    ));
                }
                let kf = permeability.permeability(1.0 - e_i);
                if !kf.is_finite() || kf <= 0.0 {
                    return Err(format!(
                        "permeability law returned {kf} at J = {:.4}",
                        1.0 - e_i
                    ));
                }
                d_node[i] = prime * kf;
            }
            // CFL: dx² / (2 · max D).
            let max_d = d_node.iter().cloned().fold(0.0f64, f64::max);
            let dt_cfl = 0.4 * dx * dx / (2.0 * max_d);
            let dt = dt_cfl.min(time - t);

            // Explicit heat-like update on the excess pressure.
            let mut dp = vec![0.0f64; n + 1];
            for i in 1..n {
                dp[i] = d_node[i] * (p[i + 1] - 2.0 * p[i] + p[i - 1]) / (dx * dx);
            }
            // Sealed platen: mirrored.
            dp[n] = d_node[n] * (p[n - 1] - p[n]) / (dx * dx);
            for i in 1..=n {
                p[i] += dt * dp[i];
                if !p[i].is_finite() {
                    return Err("excess pressure diverged".into());
                }
            }
            p[0] = 0.0;
            t += dt;
        }

        // Strain from the converged excess-pressure profile.
        (1..=n)
            .map(|i| invert_sigma_eff(sigma_eff, sigma0 - p[i]))
            .collect()
    }

    /// The creep displacement fraction `u/h = mean(e)` after marching to
    /// `time` under `σ0` (zero undrained, `σ_eff⁻¹(σ0)` at equilibrium).
    pub fn creep_with_laws(
        &self,
        sigma0: f64,
        time: f64,
        permeability: &impl PermeabilityLaw,
        sigma_eff: &impl Fn(f64) -> f64,
    ) -> Result<f64, String> {
        let e = self.strain_profile(sigma0, time, permeability, sigma_eff)?;
        Ok(e.iter().sum::<f64>() / e.len() as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creep_limits_match_analytic_bounds() {
        let m = BiphasicMaterial::default();
        let sigma = 0.2;
        // t = 0: negligible displacement (fluid carries the load). The
        // 32-term truncated series leaves a ~2% floor at t = 0.
        let initial = m.initial_displacement_fraction(sigma);
        assert!(
            initial < 0.02 * m.equilibrium_strain(sigma),
            "u(0)/u(∞) = {}",
            initial / m.equilibrium_strain(sigma)
        );
        // t → ∞: strain = σ/H_A (solid carries the load).
        let late = m.creep_displacement_fraction(sigma, 1e7, 16);
        assert!((late - m.equilibrium_strain(sigma)).abs() < 1e-9);
    }

    #[test]
    fn creep_is_monotone_and_near_gel_time_halfway() {
        let m = BiphasicMaterial::default();
        let sigma = 0.1;
        let tg = m.gel_time();
        let at_gel = m.creep_displacement_fraction(sigma, tg, 32);
        let equilibrium = m.equilibrium_strain(sigma);
        // Monotone in time.
        let quarter = m.creep_displacement_fraction(sigma, tg * 0.25, 32);
        let half = m.creep_displacement_fraction(sigma, tg * 0.5, 32);
        assert!(quarter < half && half < at_gel && at_gel < equilibrium);
        // At t = t½ the response is halfway to equilibrium (within the
        // multi-mode correction).
        assert!(
            (at_gel / equilibrium - 0.5).abs() < 0.12,
            "u(t½)/u(∞) = {}",
            at_gel / equilibrium
        );
    }

    #[test]
    fn gel_time_scales_with_thickness_squared() {
        let m = BiphasicMaterial::default();
        let doubled = BiphasicMaterial {
            thickness: m.thickness * 2.0,
            ..m
        };
        assert!((doubled.gel_time() - 4.0 * m.gel_time()).abs() < 1e-9);
    }

    #[test]
    fn shear_is_the_solid_matrix_alone() {
        let m = BiphasicMaterial::default();
        // ν_s = 0: G = H_A/2, exactly, and independent of permeability —
        // the fluid never engages in shear.
        assert!((m.solid_shear_modulus().expect("valid ν") - 0.35).abs() < 1e-12);
        let stiffer = BiphasicMaterial {
            permeability: m.permeability * 100.0,
            ..m
        };
        assert_eq!(
            m.solid_shear_modulus(),
            stiffer.solid_shear_modulus(),
            "permeability must not enter the shear response"
        );
        // General ν_s: G = H_A(1−2ν)/(2(1−ν)), hand-checked at ν = 0.3.
        let nu = BiphasicMaterial {
            poissons_ratio: 0.3,
            ..m
        };
        let expected = 0.7 * (1.0 - 0.6) / (2.0 * 0.7);
        assert!((nu.solid_shear_modulus().expect("valid ν") - expected).abs() < 1e-12);
        // ν_s ≥ ½: the linear theory has no shear modulus.
        assert!(BiphasicMaterial {
            poissons_ratio: 0.5,
            ..m
        }
        .solid_shear_modulus()
        .is_none());
    }

    #[test]
    fn unconfined_equilibrium_matches_the_solid_matrix_modulus() {
        let m = BiphasicMaterial::default();
        // ν_s = 0: the unconfined equilibrium modulus is exactly H_A —
        // the confined and unconfined long-time responses coincide.
        assert!((m.unconfined_equilibrium_modulus().expect("valid ν") - 0.7).abs() < 1e-12);
        // General ν_s: hand-checked E_s at ν = 0.25:
        // H_A · 1.25 · 0.5 / 0.75.
        let quarter = BiphasicMaterial {
            poissons_ratio: 0.25,
            ..m
        };
        let expected = 0.7 * 1.25 * 0.5 / 0.75;
        assert!(
            (quarter.unconfined_equilibrium_modulus().expect("valid ν") - expected).abs() < 1e-12
        );
        assert!(BiphasicMaterial {
            poissons_ratio: 0.5,
            ..quarter
        }
        .unconfined_equilibrium_modulus()
        .is_none());
        // For ν_s > 0 the unconfined modulus is softer than H_A, so the
        // equilibrium strain exceeds the confined one (at ν_s = 0 the two
        // coincide, as asserted above).
        let confined_strain = 0.2 / 0.7;
        let unconfined_strain = 0.2 / quarter.unconfined_equilibrium_modulus().expect("valid ν");
        assert!(unconfined_strain > confined_strain);
    }

    #[test]
    fn strain_dependent_permeability_evaluates_and_validates() {
        let material = BiphasicMaterial::default();
        // Exponential compaction law, caller-cited: k falls as the matrix
        // compacts (J < 1).
        let law = |j: f64| 0.002 * (-2.0 * (1.0 - j)).exp();
        let sweep = StrainDependentPermeability { material, law };
        // Unstrained: the law at J = 1 is the reference k.
        assert!((sweep.permeability_at(1.0).expect("valid") - 0.002).abs() < 1e-12);
        // Compaction decreases k monotonically.
        let at_half = sweep.permeability_at(0.5).expect("valid");
        assert!(at_half < 0.002 && at_half > 0.0);
        assert!(sweep.permeability_at(0.4).expect("valid") < at_half);
        // Equilibrium point: σ0 = 0.2 → J̄ = 1 − σ0/H_A = 1 − 0.2857.
        let eq = sweep.permeability_at_equilibrium(0.2).expect("valid");
        let expected_j = 1.0 - 0.2 / 0.7;
        assert!((eq - 0.002 * (-2.0f64 * (1.0 - expected_j)).exp()).abs() < 1e-12);
        // Non-positive or non-finite returns are rejected.
        assert!(StrainDependentPermeability {
            material,
            law: |_: f64| 0.0
        }
        .permeability_at(0.9)
        .is_err());
        // The constant baseline reproduces the material's own k at any J.
        let constant = StrainDependentPermeability {
            material,
            law: ConstantPermeability { k: 0.002 },
        };
        assert_eq!(
            constant.permeability_at(0.7).expect("valid"),
            constant.material.permeability
        );
    }

    #[test]
    fn unconfined_instantaneous_modulus_brackets_the_equilibrium() {
        let m = BiphasicMaterial::default();
        // ν_s = 0: E(0⁺) = 3G = 1.5·H_A, against the E(∞) = H_A bookend.
        assert!((m.unconfined_instantaneous_modulus().expect("valid ν") - 1.05).abs() < 1e-12);
        // General ν_s: the 3G identity, hand-checked at ν = 0.3:
        // 3 · H_A(1−2ν)/(2(1−ν)) = 6H_A/7.
        let nu = BiphasicMaterial {
            poissons_ratio: 0.3,
            ..m
        };
        assert!(
            (nu.unconfined_instantaneous_modulus().expect("valid ν") - 6.0 * 0.7 / 7.0).abs()
                < 1e-12
        );
        // The relaxation bookends are ordered: E(0⁺) > E(∞) for ν_s > 0,
        // in the ratio 3/(2(1+ν_s)) — at ν_s = 0 the instantaneous wall
        // is still 1.5× the equilibrium one.
        let e_inst = nu.unconfined_instantaneous_modulus().expect("valid ν");
        let e_eq = nu.unconfined_equilibrium_modulus().expect("valid ν");
        assert!((e_inst / e_eq - 3.0 / (2.0 * 1.3)).abs() < 1e-12);
        assert!(e_inst > e_eq);
        assert!(
            (m.unconfined_instantaneous_modulus().expect("valid ν")
                / m.unconfined_equilibrium_modulus().expect("valid ν")
                - 1.5)
                .abs()
                < 1e-12
        );
        // ν_s ≥ ½: no linear instantaneous modulus.
        assert!(BiphasicMaterial {
            poissons_ratio: 0.5,
            ..m
        }
        .unconfined_instantaneous_modulus()
        .is_none());
    }

    #[test]
    fn stefan_closed_form_matches_its_defining_ode() {
        let film = SqueezeFilm::new(3.0e-8, 12.0).expect("valid");
        // 5 N keeps the RK4 reference inside its stability region at this
        // step size (the film's rate constant scales with the load).
        let (h0, load) = (0.5, 5.0);
        // Integrate dh/dt = −2Wh³/(3πμR⁴) (Stefan's load relation solved
        // for the rate) numerically, and compare against the closed form.
        let rate = |h: f64| {
            -2.0 * load * h.powi(3)
                / (3.0 * core::f64::consts::PI * film.viscosity * film.radius.powi(4))
        };
        let mut h = h0;
        let t_end = 50.0;
        let steps = 200_000;
        let dt = t_end / steps as f64;
        for _ in 0..steps {
            // RK4 on the scalar ODE.
            let k1 = rate(h);
            let k2 = rate(h + 0.5 * dt * k1);
            let k3 = rate(h + 0.5 * dt * k2);
            let k4 = rate(h + dt * k3);
            h += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }
        let closed = film.film_thickness(h0, load, t_end);
        assert!(
            (h - closed).abs() < 1e-4 * closed,
            "RK4 {h} vs closed form {closed}"
        );
        // The exact squeeze time inverts the closed form.
        let t_half = film.time_to_squeeze(load, h0, closed).expect("thinning");
        assert!(
            (film.film_thickness(h0, load, t_half) - closed).abs() < 1e-9 * closed,
            "time inverse must land on the same thickness"
        );
    }

    #[test]
    fn squeeze_film_limits_and_load_capacity_round_trip() {
        let film = SqueezeFilm::new(3.0e-8, 12.0).expect("valid");
        let (h0, load) = (0.5, 500.0);
        // h(0) = h0 exactly; the film thins monotonically and never
        // touches; more load or more time squeezes further.
        assert!((film.film_thickness(h0, load, 0.0) - h0).abs() < 1e-12);
        let t1 = film.film_thickness(h0, load, 10.0);
        let t2 = film.film_thickness(h0, load, 20.0);
        assert!(h0 > t1 && t1 > t2 && t2 > 0.0);
        assert!(film.film_thickness(h0, 2.0 * load, 10.0) < t1);
        // Load capacity round-trips: the approach rate the closed form
        // implies at thickness h carries exactly the applied load.
        let h = t1;
        let dh_dt = (film.film_thickness(h0, load, 10.0 + 1e-6)
            - film.film_thickness(h0, load, 10.0))
            / 1e-6;
        let capacity = film.load_capacity(h, -dh_dt);
        assert!(
            (capacity - load).abs() < 1e-3 * load,
            "capacity {capacity} vs load {load}"
        );
        // Construction validation.
        assert!(SqueezeFilm::new(0.0, 12.0).is_err());
        assert!(SqueezeFilm::new(-1.0, 12.0).is_err());
        assert!(SqueezeFilm::new(3.0e-8, 0.0).is_err());
        assert!(SqueezeFilm::new(f64::NAN, 12.0).is_err());
        // Time-to-squeeze rejects non-thinning and non-positive inputs.
        assert!(film.time_to_squeeze(load, 0.5, 0.5).is_err());
        assert!(film.time_to_squeeze(load, 0.5, 0.6).is_err());
        assert!(film.time_to_squeeze(load, 0.0, 0.1).is_err());
        assert!(film.time_to_squeeze(-1.0, 0.5, 0.1).is_err());
    }

    #[test]
    fn constant_k_stepper_reproduces_the_series_solution() {
        let m = BiphasicMaterial::default();
        let stepper = ConfinedCreepStepper::new(m, 150).expect("valid");
        let law = ConstantPermeability { k: m.permeability };
        let sigma = 0.1;
        let equilibrium = m.equilibrium_strain(sigma);
        // t = 0: the fluid carries everything, the displacement is zero.
        let u0 = stepper.creep_fraction_with_law(sigma, 1e-6, &law);
        assert!(
            u0.abs() < 0.01 * equilibrium,
            "u(0)/u(∞) = {}",
            u0 / equilibrium
        );
        // Mid-transient: the numerics reproduce the closed form to a
        // small fraction of the equilibrium strain (2nd-order space,
        // 1st-order time, 200 cells).
        let tg = m.gel_time();
        for frac in [0.25, 0.5, 1.0, 2.0] {
            let numeric = stepper.creep_fraction_with_law(sigma, tg * frac, &law);
            let series = m.creep_displacement_fraction(sigma, tg * frac, 32);
            assert!(
                (numeric - series).abs() < 0.01 * equilibrium,
                "t = {frac}·t½: numeric {numeric} vs series {series}"
            );
        }
        // Long time: equilibrium.
        let late = stepper.creep_fraction_with_law(sigma, tg * 20.0, &law);
        assert!((late - equilibrium).abs() < 1e-3 * equilibrium);
    }

    #[test]
    fn stepper_error_shrinks_under_refinement() {
        let m = BiphasicMaterial::default();
        let law = ConstantPermeability { k: m.permeability };
        let sigma = 0.1;
        let t = m.gel_time() * 0.5;
        let series = m.creep_displacement_fraction(sigma, t, 32);
        let error = |cells: usize| {
            let s = ConfinedCreepStepper::new(m, cells).expect("valid");
            (s.creep_fraction_with_law(sigma, t, &law) - series).abs()
        };
        let (e100, e200) = (error(100), error(200));
        assert!(
            e200 < e100 * 0.7,
            "refinement must shrink the error: {e100} -> {e200}"
        );
    }

    #[test]
    fn strain_dependent_permeability_slows_the_creep() {
        let m = BiphasicMaterial::default();
        let stepper = ConfinedCreepStepper::new(m, 100).expect("valid");
        let constant = ConstantPermeability { k: m.permeability };
        // The compaction law evaluated at J = 1 equals the reference k, so
        // the two runs start identically; as the matrix compacts the law's
        // k falls, drainage slows, and the strain at any finite time is
        // smaller than the constant-k run.
        let law = |j: f64| m.permeability * (-2.0 * (1.0 - j)).exp();
        let sigma = 0.1;
        let t = m.gel_time();
        let u_constant = stepper.creep_fraction_with_law(sigma, t, &constant);
        let u_compacting = stepper.creep_fraction_with_law(sigma, t, &law);
        assert!(
            u_compacting < u_constant,
            "compacted matrix must drain slower: {u_compacting} vs {u_constant}"
        );
        // Both reach the same equilibrium eventually (the compacting run
        // needs longer: its drained-state k is ~1.3x smaller).
        let late_c = stepper.creep_fraction_with_law(sigma, t * 30.0, &constant);
        let late_l = stepper.creep_fraction_with_law(sigma, t * 30.0, &law);
        let equilibrium = m.equilibrium_strain(sigma);
        assert!((late_c - equilibrium).abs() < 1e-3 * equilibrium);
        assert!((late_l - equilibrium).abs() < 1e-3 * equilibrium);
    }

    #[test]
    fn stepper_validates_its_construction() {
        let m = BiphasicMaterial::default();
        assert!(ConfinedCreepStepper::new(m, 4).is_err());
        assert!(ConfinedCreepStepper::new(m, 8).is_ok());
    }

    #[test]
    fn fluid_pressure_decays() {
        let m = BiphasicMaterial::default();
        let early = m.fluid_pressure_fraction(0.05, 64);
        let late = m.fluid_pressure_fraction(m.gel_time() * 10.0, 64);
        assert!(early > 0.1, "early pressure {early}");
        assert!(late < early);
        assert!(late < 0.05, "late pressure {late}");
    }

    // ---- Nonlinear solid-matrix stepper ----

    #[test]
    fn nonlinear_linear_law_reproduces_the_linear_stepper() {
        // The cross-verification: sigma_eff = H_A (J - 1) makes the
        // nonlinear stepper's diffusivity k*H_A and its Dirichlet stretch
        // 1 + sigma0/H_A — the linear ConfinedCreepStepper's exact
        // problem. The two independent steppers must agree.
        let m = BiphasicMaterial::default();
        let lin = ConfinedCreepStepper::new(m, 100).expect("valid");
        let nl = NonlinearConfinedStepper::new(m.thickness, 100).expect("valid");
        let sigma0 = 0.1;
        let linear_law = ConstantPermeability { k: m.permeability };
        // In strain form: sigma_eff(e) = H_A e — exactly the linear law.
        let sigma_eff = |e: f64| m.aggregate_modulus * e;
        let tg = m.gel_time();
        for frac in [0.5, 1.0, 2.0] {
            let from_linear = lin.creep_fraction_with_law(sigma0, tg * frac, &linear_law);
            let from_nonlinear = nl
                .creep_with_laws(sigma0, tg * frac, &linear_law, &sigma_eff)
                .expect("marches");
            assert!(
                (from_linear - from_nonlinear).abs() < 0.01 * sigma0,
                "t = {frac} t1/2: linear stepper {from_linear} vs nonlinear {from_nonlinear}"
            );
        }
    }

    #[test]
    fn nonlinear_equilibrium_matches_the_inverted_law() {
        let m = BiphasicMaterial::default();
        let nl = NonlinearConfinedStepper::new(m.thickness, 100).expect("valid");
        let linear_law = ConstantPermeability { k: m.permeability };
        let sigma0 = 0.1;
        // An exponentially stiffening solid matrix (caller-cited form):
        // sigma_eff(e) = 0.3 (e^{0.5 e} - 1); its equilibrium strain is
        // sigma_eff^{-1}(sigma0), computed here by the same bisection the
        // stepper uses.
        let m_stiff = 0.5f64;
        let sigma_eff = |e: f64| (m_stiff * e).exp_m1() * 0.3;
        let j_eq = {
            let (mut lo, mut hi) = (0.05, 2.0);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if sigma_eff(mid) < sigma0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        };
        let late = nl
            .creep_with_laws(sigma0, 3.0e4, &linear_law, &sigma_eff)
            .expect("marches");
        let expected = j_eq;
        assert!(
            (late - expected).abs() < 1e-3 * expected.abs().max(1e-6),
            "late {late} vs equilibrium strain {expected}"
        );
    }

    #[test]
    fn nonlinear_stiffening_law_slows_the_creep() {
        let m = BiphasicMaterial::default();
        let nl = NonlinearConfinedStepper::new(m.thickness, 100).expect("valid");
        let linear_law = ConstantPermeability { k: m.permeability };
        let sigma0 = 0.1;
        let linear = |e: f64| 0.7 * e;
        // A monotone law that stiffens in compression: dSigma/de = 0.7 +
        // 6 e^2 > 0 everywhere, and the stress rises faster than linear
        // as e grows, so the same sigma0 needs less strain at equilibrium
        // and the transient is slower.
        let stiffening = |e: f64| 0.7 * e + 2.0 * e.powi(3);
        let t = m.gel_time();
        let u_lin = nl
            .creep_with_laws(sigma0, t, &linear_law, &linear)
            .expect("marches");
        let u_stiff = nl
            .creep_with_laws(sigma0, t, &linear_law, &stiffening)
            .expect("marches");
        assert!(
            u_stiff < u_lin,
            "stiffening law must creep less at the same time: {u_stiff} vs {u_lin}"
        );
    }

    #[test]
    fn nonlinear_stepper_validates_and_inverts() {
        let nl = NonlinearConfinedStepper::new(2.0, 4);
        assert!(nl.is_err());
        let nl = NonlinearConfinedStepper::new(2.0, 8).expect("valid");
        let linear_law = ConstantPermeability { k: 0.002 };
        let non_monotone = |e: f64| 0.7 * e * e; // U-shaped: not increasing
        assert!(
            nl.creep_with_laws(0.1, 1.0, &linear_law, &non_monotone)
                .is_err(),
            "a non-monotone law must be rejected, not silently marched"
        );
    }
}
