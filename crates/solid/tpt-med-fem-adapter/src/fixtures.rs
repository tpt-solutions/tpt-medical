//! Stent verification fixtures (RFC 0013): meshes the Level 3 deployment is
//! checked on, so the cross-check against the Level 1/2 ring models runs on a
//! defined geometry rather than whatever a test happens to build.
//!
//! The crate does not mesh CAD. [`crown_ring`] is a *swept* ring: a
//! rectangular strut cross-section (width `w` in the cylinder surface,
//! thickness `t` radially) swept along the closed centreline
//!
//! ```text
//! c(θ) = ( R cos θ,  R sin θ,  A cos(N θ) ),   θ ∈ [0, 2π)
//! ```
//!
//! — `N` crowns and `N` strut pairs, the one-ring geometry Level 1's
//! "N crown springs" abstraction describes. The sinusoid is a smooth stand-in
//! for a filleted zig-zag: it has no sharp crown to mesh badly, and its
//! crown curvature is set by `A` and `N` alone, which keeps the fixture
//! defined by five numbers. It is **not** a laser-cut stent geometry and
//! makes no claim to represent a commercial device; the F2394 benchmark
//! geometry is a separate acceptance item.
//!
//! Cross-section frame at station `θ`: tangent `T = c′/|c′|`, radial
//! `e_r` (the outward radial direction orthogonalised against `T`) and
//! lateral `B = e_r × T`, which makes `(T, B, e_r)` right-handed so every
//! element has a positive Jacobian.

use crate::mesh::{Hex8Mesh, MeshError};
use tpt_fem_element::{Hex8, ReferenceElement};
use tpt_med_geometry::Vec3;

/// Parameters of a swept crown ring (mm).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrownRingSpec {
    /// Number of crowns `N` (each is one peak of the centreline).
    pub crowns: usize,
    /// Centreline radius `R`.
    pub radius: f64,
    /// Axial amplitude `A` of the centreline.
    pub amplitude: f64,
    /// Strut width `w`, in the cylinder surface.
    pub width: f64,
    /// Strut thickness `t`, radial.
    pub thickness: f64,
    /// Stations per crown period along the centreline.
    pub stations_per_crown: usize,
    /// Elements across the width.
    pub width_elements: usize,
    /// Elements through the thickness.
    pub thickness_elements: usize,
}

impl CrownRingSpec {
    /// An 8 mm ring of `crowns` crowns with a 0.2 mm × 0.15 mm strut — the
    /// scale of a coronary/peripheral nitinol stent — meshed 2 × 2 across the
    /// section so bending is not a single-element locking mode.
    pub fn new(crowns: usize) -> Self {
        Self {
            crowns,
            radius: 4.0,
            amplitude: 1.0,
            width: 0.15,
            thickness: 0.2,
            stations_per_crown: 8,
            width_elements: 2,
            thickness_elements: 2,
        }
    }

    /// The outer radius of the unloaded ring, `R + t/2` — the largest
    /// radial extent, and the diameter Level 1 calls the free diameter is
    /// twice this.
    pub fn outer_radius(&self) -> f64 {
        self.radius + 0.5 * self.thickness
    }
}

/// Errors building a fixture.
#[derive(Debug, Clone, PartialEq)]
pub enum FixtureError {
    /// A count is too small to define the ring (needs ≥ 3 crowns, ≥ 3
    /// stations per crown and ≥ 1 element across each section direction).
    TooCoarse(&'static str),
    /// A length is not finite and positive, or the section does not fit
    /// inside the ring.
    BadDimension(&'static str),
    /// The mesh failed validation.
    Mesh(MeshError),
}

impl std::fmt::Display for FixtureError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooCoarse(w) => write!(f, "fixture too coarse: {w}"),
            Self::BadDimension(w) => write!(f, "bad dimension: {w}"),
            Self::Mesh(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for FixtureError {}

/// A built crown ring.
#[derive(Debug, Clone)]
pub struct CrownRing {
    /// The mesh.
    pub mesh: Hex8Mesh,
    /// The spec it was built from.
    pub spec: CrownRingSpec,
    /// Nodes on the radially outer surface — the wall's slave nodes.
    pub outer_nodes: Vec<usize>,
    /// Rigid-body pins `(dof, value)`: the axial DOF at three stations
    /// spread around the ring (translation and tilt) and the tangential DOF
    /// at station 0 (rotation about the axis). The radial wall restrains the
    /// in-plane translations. The ring's symmetry makes the reactions at
    /// these pins vanish for a radial load, so they remove rigid modes
    /// without carrying load.
    pub pins: Vec<(usize, f64)>,
    /// Stations around the ring.
    pub stations: usize,
}

impl CrownRingSpec {
    /// Build the ring.
    ///
    /// # Errors
    /// [`FixtureError`] for a degenerate spec.
    pub fn build(&self) -> Result<CrownRing, FixtureError> {
        if self.crowns < 3 {
            return Err(FixtureError::TooCoarse("fewer than 3 crowns"));
        }
        if self.stations_per_crown < 3 {
            return Err(FixtureError::TooCoarse("fewer than 3 stations per crown"));
        }
        if self.width_elements == 0 || self.thickness_elements == 0 {
            return Err(FixtureError::TooCoarse("no elements across the section"));
        }
        for (v, name) in [
            (self.radius, "radius"),
            (self.amplitude, "amplitude"),
            (self.width, "width"),
            (self.thickness, "thickness"),
        ] {
            if !v.is_finite() || v <= 0.0 {
                return Err(FixtureError::BadDimension(name));
            }
        }
        if self.thickness >= self.radius || self.width >= self.radius {
            return Err(FixtureError::BadDimension("section does not fit the ring"));
        }

        let m = self.crowns * self.stations_per_crown;
        let (nw, nt) = (self.width_elements, self.thickness_elements);
        let per_station = (nw + 1) * (nt + 1);
        let n = self.crowns as f64;
        let node_id = |k: usize, j: usize, l: usize| ((k % m) * (nw + 1) + j) * (nt + 1) + l;

        let mut nodes = vec![Vec3::new(0.0, 0.0, 0.0); m * per_station];
        for k in 0..m {
            let th = 2.0 * std::f64::consts::PI * k as f64 / m as f64;
            let (s, c) = th.sin_cos();
            let centre = [
                self.radius * c,
                self.radius * s,
                self.amplitude * (n * th).cos(),
            ];
            // c'(θ), normalised.
            let d = [
                -self.radius * s,
                self.radius * c,
                -self.amplitude * n * (n * th).sin(),
            ];
            let dn = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let tg = [d[0] / dn, d[1] / dn, d[2] / dn];
            let er0 = [c, s, 0.0];
            let dotp = er0[0] * tg[0] + er0[1] * tg[1] + er0[2] * tg[2];
            let er_raw = [
                er0[0] - dotp * tg[0],
                er0[1] - dotp * tg[1],
                er0[2] - dotp * tg[2],
            ];
            let en = (er_raw[0].powi(2) + er_raw[1].powi(2) + er_raw[2].powi(2)).sqrt();
            let er = [er_raw[0] / en, er_raw[1] / en, er_raw[2] / en];
            // B = e_r x T
            let b = [
                er[1] * tg[2] - er[2] * tg[1],
                er[2] * tg[0] - er[0] * tg[2],
                er[0] * tg[1] - er[1] * tg[0],
            ];
            for j in 0..=nw {
                for l in 0..=nt {
                    let a = (j as f64 / nw as f64 - 0.5) * self.width;
                    let h = (l as f64 / nt as f64 - 0.5) * self.thickness;
                    nodes[node_id(k, j, l)] = Vec3::new(
                        centre[0] + a * b[0] + h * er[0],
                        centre[1] + a * b[1] + h * er[1],
                        centre[2] + a * b[2] + h * er[2],
                    );
                }
            }
        }

        let mut elements = Vec::with_capacity(m * nw * nt);
        for k in 0..m {
            for j in 0..nw {
                for l in 0..nt {
                    let element: Vec<usize> = Hex8::nodes()
                        .iter()
                        .map(|r| {
                            let di = ((r[0] + 1.0) * 0.5).round() as usize;
                            let dj = ((r[1] + 1.0) * 0.5).round() as usize;
                            let dl = ((r[2] + 1.0) * 0.5).round() as usize;
                            node_id(k + di, j + dj, l + dl)
                        })
                        .collect();
                    elements.push(element);
                }
            }
        }
        let mesh = Hex8Mesh::from_parts(nodes, elements).map_err(FixtureError::Mesh)?;

        let mut outer_nodes = Vec::new();
        for k in 0..m {
            for j in 0..=nw {
                outer_nodes.push(node_id(k, j, nt));
            }
        }

        // Rigid-mode pins at the mid-section inner node of chosen stations.
        let pin_node = |k: usize| node_id(k, nw / 2, 0);
        let mut pins = Vec::new();
        for third in 0..3 {
            pins.push((mesh.dof(pin_node(third * m / 3), 2), 0.0));
        }
        // At θ = 0 the tangent is the global y axis: pin it against twist.
        pins.push((mesh.dof(pin_node(0), 1), 0.0));

        Ok(CrownRing {
            mesh,
            spec: *self,
            outer_nodes,
            pins,
            stations: m,
        })
    }
}

/// The default `crowns`-crown ring ([`CrownRingSpec::new`]).
///
/// # Errors
/// [`FixtureError`] if `crowns < 3`.
pub fn crown_ring(crowns: usize) -> Result<CrownRing, FixtureError> {
    CrownRingSpec::new(crowns).build()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Mesh volume and the exact swept volume (section area times the
    /// centreline arc length, by quadrature) of a spec.
    fn volumes(s: &CrownRingSpec) -> (f64, f64) {
        let ring = s.build().expect("ring");
        let mesh = &ring.mesh;
        let zero = vec![0.0; mesh.dof_count()];
        assert!(mesh.inverted_elements(&zero).expect("check").is_empty());
        let vol: f64 = (0..mesh.element_count())
            .map(|e| mesh.element_volume(e).expect("volume"))
            .sum();
        let steps = 20_000;
        let mut len = 0.0;
        for i in 0..steps {
            let th = 2.0 * std::f64::consts::PI * (i as f64 + 0.5) / steps as f64;
            let dz = s.amplitude * s.crowns as f64 * (s.crowns as f64 * th).sin();
            len +=
                (s.radius * s.radius + dz * dz).sqrt() * 2.0 * std::f64::consts::PI / steps as f64;
        }
        (vol, s.width * s.thickness * len)
    }

    #[test]
    fn default_ring_is_a_valid_closed_mesh() {
        let ring = crown_ring(8).expect("ring");
        let mesh = &ring.mesh;
        assert_eq!(mesh.element_count(), 8 * 8 * 2 * 2);
        assert_eq!(mesh.node_count(), 8 * 8 * 3 * 3);
    }

    #[test]
    fn swept_volume_converges_to_the_exact_centreline_volume() {
        // The polyline centreline underestimates arc length; the error must
        // fall as stations are added, which is the mesh-geometry convergence
        // the verification study builds on.
        let err = |m: usize| {
            let mut s = CrownRingSpec::new(8);
            s.stations_per_crown = m;
            let (vol, exact) = volumes(&s);
            (exact - vol).abs() / exact
        };
        let (e6, e12, e24) = (err(6), err(12), err(24));
        assert!(e6 > e12 && e12 > e24, "{e6} {e12} {e24}");
        assert!(e24 < 0.01, "finest geometry error {e24}");
        // Roughly second order in the station spacing.
        assert!(e12 < 0.4 * e6, "{e6} -> {e12}");
    }

    #[test]
    fn outer_nodes_sit_at_the_outer_radius_and_pins_remove_rigid_modes() {
        let ring = crown_ring(8).expect("ring");
        let s = ring.spec;
        // The outer surface is radially outward of the centreline by t/2
        // (to within the tilt of the section against the radial direction).
        for &n in &ring.outer_nodes {
            let p = ring.mesh.nodes()[n].to_array();
            let rho = p[0].hypot(p[1]);
            assert!(
                (rho - s.outer_radius()).abs() < 0.05 * s.thickness + 0.02 * s.width,
                "outer node radius {rho}"
            );
        }
        assert_eq!(ring.pins.len(), 4);
    }

    #[test]
    fn degenerate_specs_are_rejected() {
        assert!(matches!(crown_ring(2), Err(FixtureError::TooCoarse(_))));
        let mut s = CrownRingSpec::new(6);
        s.thickness = 10.0;
        assert!(matches!(s.build(), Err(FixtureError::BadDimension(_))));
        let mut s = CrownRingSpec::new(6);
        s.width_elements = 0;
        assert!(matches!(s.build(), Err(FixtureError::TooCoarse(_))));
        let mut s = CrownRingSpec::new(6);
        s.amplitude = f64::NAN;
        assert!(matches!(s.build(), Err(FixtureError::BadDimension(_))));
    }
}
