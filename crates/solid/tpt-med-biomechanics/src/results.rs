//! Solver results and scalar stress measures.

use tpt_med_geometry::Vec3;

/// Per-element stress state (Voigt `[xx, yy, zz, xy, yz, xz]`, MPa).
#[derive(Debug, Clone, Copy)]
pub struct ElementStress {
    /// Element id.
    pub element: u32,
    /// Strain components (engineering shear).
    pub strain: [f64; 6],
    /// Stress components (MPa).
    pub stress: [f64; 6],
}

impl ElementStress {
    /// von Mises equivalent stress (MPa).
    pub fn von_mises(&self) -> f64 {
        let s = &self.stress;
        let (sx, sy, sz) = (s[0], s[1], s[2]);
        let (txy, tyz, txz) = (s[3], s[4], s[5]);
        (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2))
            + 3.0 * (txy * txy + tyz * tyz + txz * txz))
            .sqrt()
    }

    /// Principal stresses sorted descending (analytic cubic, MPa).
    pub fn principal_stresses(&self) -> [f64; 3] {
        let s = self.stress;
        let m = tpt_med_geometry::Mat3::from_array([
            s[0], s[3], s[5], s[3], s[1], s[4], s[5], s[4], s[2],
        ]);
        m.symmetric_eigenvalues()
    }

    /// Hydrostatic (mean) stress, MPa.
    pub fn hydrostatic(&self) -> f64 {
        (self.stress[0] + self.stress[1] + self.stress[2]) / 3.0
    }
}

/// Complete static solution.
#[derive(Debug, Clone)]
pub struct StressResult {
    /// Nodal displacement field (mm).
    pub displacements: Vec<Vec3>,
    /// Linear-solver statistics.
    pub stats: crate::sparse::SolveStats,
    /// Per-element stress states (MPa).
    pub stresses: Vec<ElementStress>,
}

impl StressResult {
    /// Maximum von Mises stress over all elements (MPa).
    pub fn max_von_mises(&self) -> f64 {
        self.stresses
            .iter()
            .map(|e| e.von_mises())
            .fold(0.0, f64::max)
    }

    /// Mean von Mises stress (MPa).
    pub fn mean_von_mises(&self) -> f64 {
        if self.stresses.is_empty() {
            return 0.0;
        }
        self.stresses.iter().map(|e| e.von_mises()).sum::<f64>() / self.stresses.len() as f64
    }

    /// Maximum displacement magnitude (mm).
    pub fn max_displacement(&self) -> f64 {
        self.displacements
            .iter()
            .map(|d| d.norm())
            .fold(0.0, f64::max)
    }

    /// Element id with the highest von Mises stress.
    pub fn critical_element(&self) -> Option<u32> {
        self.stresses
            .iter()
            .max_by(|a, b| {
                a.von_mises()
                    .partial_cmp(&b.von_mises())
                    .unwrap_or(core::cmp::Ordering::Equal)
            })
            .map(|e| e.element)
    }
}
