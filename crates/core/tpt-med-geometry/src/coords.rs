//! Anatomical and imaging coordinate systems.
//!
//! DICOM patient coordinates use **LPS** (x: Left, y: Posterior, z: Superior)
//! while research tools (NIfTI, 3D Slicer) use **RAS** (x: Right, y:
//! Anterior, z: Superior). Conversions flip the x and y axes. Anatomical
//! direction conventions follow the International Society of Biomechanics
//! (ISB) recommendations (Wu & Cavanagh 1995; Wu et al. 2002).

use crate::mat3::Mat3;
use crate::vec3::Vec3;

/// A named coordinate system relevant to the imaging/simulation pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateSystem {
    /// DICOM patient coordinate system (+x Left, +y Posterior, +z Superior).
    Lps,
    /// Research coordinate system (+x Right, +y Anterior, +z Superior).
    Ras,
}

/// An image frame: origin (first voxel center, patient coords), column
/// direction cosine, and row direction cosine, exactly as DICOM tags
/// (0020,0032) and (0020,0037) define them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ImageFrame {
    /// Position of the first transmitted voxel (ImagePositionPatient).
    pub origin: Vec3,
    /// Direction cosine of the image rows (first 3 values of
    /// ImageOrientationPatient).
    pub row_dir: Vec3,
    /// Direction cosine of the image columns (last 3 values of
    /// ImageOrientationPatient).
    pub col_dir: Vec3,
}

impl ImageFrame {
    /// Slice normal (row × column), pointing in the direction of increasing
    /// slice position per DICOM PS3.3 C.7.6.2.1.1.
    pub fn normal(&self) -> Vec3 {
        self.row_dir.cross(self.col_dir).normalize()
    }

    /// Patient-space position of the voxel at (col, row) in this slice.
    pub fn voxel_position(&self, col: usize, row: usize, spacing: (f64, f64)) -> Vec3 {
        self.origin
            + self.row_dir * (col as f64 * spacing.0)
            + self.col_dir * (row as f64 * spacing.1)
    }

    /// The 3×3 rotation from image indices to patient coordinates.
    pub fn rotation(&self) -> Mat3 {
        Mat3::from_cols(
            self.row_dir.normalize(),
            self.col_dir.normalize(),
            self.normal(),
        )
    }
}

/// Flips x and y: LPS → RAS or RAS → LPS (the map is an involution).
pub fn flip_xy(p: Vec3) -> Vec3 {
    Vec3::new(-p.x, -p.y, p.z)
}

/// Converts a point from DICOM LPS to research RAS.
pub fn lps_to_ras(p: Vec3) -> Vec3 {
    flip_xy(p)
}

/// Converts a point from research RAS to DICOM LPS.
pub fn ras_to_lps(p: Vec3) -> Vec3 {
    flip_xy(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lps_ras_roundtrip() {
        let p = Vec3::new(-12.5, 30.0, 105.0); // left, posterior, superior
        let ras = lps_to_ras(p);
        assert_eq!(ras, Vec3::new(12.5, -30.0, 105.0)); // right, anterior, superior
        assert_eq!(ras_to_lps(ras), p);
    }

    #[test]
    fn image_frame_normal_and_positions() {
        // Axial slice: rows along +x (LPS left→right is +? ) — standard
        // axial CT: row_dir = [1,0,0], col_dir = [0,1,0], normal = +z.
        let frame = ImageFrame {
            origin: Vec3::new(-100.0, -100.0, 5.0),
            row_dir: Vec3::X,
            col_dir: Vec3::Y,
        };
        assert_eq!(frame.normal(), Vec3::Z);
        let p = frame.voxel_position(10, 20, (0.5, 0.5));
        assert_eq!(p, Vec3::new(-95.0, -90.0, 5.0));
        // Rotation columns are orthonormal
        let r = frame.rotation();
        assert!((r.det() - 1.0).abs() < 1e-12);
    }
}
