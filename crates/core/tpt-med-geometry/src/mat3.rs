//! Row-major 3×3 matrix.

use crate::vec3::Vec3;

/// Row-major 3×3 matrix. Rows are stored as three [`Vec3`]s so indexing
/// reads naturally: `m.row(0)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat3 {
    r0: Vec3,
    r1: Vec3,
    r2: Vec3,
}

impl Mat3 {
    /// Identity matrix.
    pub const IDENTITY: Mat3 = Mat3::from_rows(
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    );

    /// Zero matrix.
    pub const ZERO: Mat3 = Mat3::from_rows(Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);

    /// Builds a matrix from three rows.
    #[inline]
    pub const fn from_rows(r0: Vec3, r1: Vec3, r2: Vec3) -> Self {
        Self { r0, r1, r2 }
    }

    /// Builds a matrix from a row-major array.
    #[inline]
    pub fn from_array(a: [f64; 9]) -> Self {
        Self::from_rows(
            Vec3::new(a[0], a[1], a[2]),
            Vec3::new(a[3], a[4], a[5]),
            Vec3::new(a[6], a[7], a[8]),
        )
    }

    /// Builds a matrix from column vectors (DICOM orientation tags supply
    /// direction cosines as columns).
    #[inline]
    pub fn from_cols(c0: Vec3, c1: Vec3, c2: Vec3) -> Self {
        Self::from_rows(
            Vec3::new(c0.x, c1.x, c2.x),
            Vec3::new(c0.y, c1.y, c2.y),
            Vec3::new(c0.z, c1.z, c2.z),
        )
    }

    /// Diagonal matrix.
    #[inline]
    pub fn diagonal(d: [f64; 3]) -> Self {
        Self::from_rows(
            Vec3::new(d[0], 0.0, 0.0),
            Vec3::new(0.0, d[1], 0.0),
            Vec3::new(0.0, 0.0, d[2]),
        )
    }

    /// Row `i`.
    #[inline]
    pub fn row(&self, i: usize) -> Vec3 {
        match i {
            0 => self.r0,
            1 => self.r1,
            _ => self.r2,
        }
    }

    /// Column `i`.
    #[inline]
    pub fn col(&self, i: usize) -> Vec3 {
        Vec3::new(
            self.r0.to_array()[i],
            self.r1.to_array()[i],
            self.r2.to_array()[i],
        )
    }

    /// Element at `(row, col)`.
    #[inline]
    pub fn at(&self, r: usize, c: usize) -> f64 {
        self.row(r).to_array()[c]
    }

    /// Sets the element at `(row, col)`.
    #[inline]
    pub fn set(&mut self, r: usize, c: usize, v: f64) {
        let arr = match r {
            0 => &mut self.r0,
            1 => &mut self.r1,
            _ => &mut self.r2,
        };
        let mut a = arr.to_array();
        a[c] = v;
        *arr = Vec3::new(a[0], a[1], a[2]);
    }

    /// Row-major flat array.
    pub fn to_array(self) -> [f64; 9] {
        [
            self.r0.x, self.r0.y, self.r0.z, self.r1.x, self.r1.y, self.r1.z, self.r2.x, self.r2.y,
            self.r2.z,
        ]
    }

    /// Matrix–vector product.
    #[inline]
    pub fn mul_vec(&self, v: Vec3) -> Vec3 {
        Vec3::new(self.r0.dot(v), self.r1.dot(v), self.r2.dot(v))
    }

    /// Matrix product `self * rhs`.
    pub fn mul_mat(&self, rhs: &Mat3) -> Mat3 {
        Mat3::from_cols(
            self.mul_vec(rhs.col(0)),
            self.mul_vec(rhs.col(1)),
            self.mul_vec(rhs.col(2)),
        )
    }

    /// Transpose.
    pub fn transpose(&self) -> Mat3 {
        Mat3::from_cols(self.r0, self.r1, self.r2)
    }

    /// Determinant.
    pub fn det(&self) -> f64 {
        self.r0.dot(self.r1.cross(self.r2))
    }

    /// Inverse; returns `None` for (near-)singular matrices.
    pub fn inverse(&self) -> Option<Mat3> {
        let d = self.det();
        if d.abs() < crate::EPS {
            return None;
        }
        // Columns of the inverse are cross products of row pairs, divided
        // by the determinant (adjugate / det).
        let c0 = self.r1.cross(self.r2) / d;
        let c1 = self.r2.cross(self.r0) / d;
        let c2 = self.r0.cross(self.r1) / d;
        Some(Mat3::from_cols(c0, c1, c2))
    }

    /// Rotation about an arbitrary axis (right-hand rule), radians.
    pub fn rotation_axis_angle(axis: Vec3, angle: f64) -> Mat3 {
        let k = axis.normalize();
        let (s, c) = angle.sin_cos();
        let c1 = 1.0 - c;
        let (kx, ky, kz) = (k.x, k.y, k.z);
        Mat3::from_array([
            c + kx * kx * c1,
            kx * ky * c1 - kz * s,
            kx * kz * c1 + ky * s,
            ky * kx * c1 + kz * s,
            c + ky * ky * c1,
            ky * kz * c1 - kx * s,
            kz * kx * c1 - ky * s,
            kz * ky * c1 + kx * s,
            c + kz * kz * c1,
        ])
    }

    /// Scaling matrix.
    pub fn scaling(s: Vec3) -> Mat3 {
        Mat3::diagonal([s.x, s.y, s.z])
    }

    /// Trace.
    pub fn trace(&self) -> f64 {
        self.at(0, 0) + self.at(1, 1) + self.at(2, 2)
    }

    /// Eigenvalues of a **symmetric** matrix, sorted descending, via the
    /// trigonometric solution of the characteristic cubic. Used for
    /// principal stresses and strains. Only meaningful for symmetric input.
    pub fn symmetric_eigenvalues(&self) -> [f64; 3] {
        // Characteristic polynomial: λ³ - I1 λ² + I2 λ - I3 = 0
        let i1 = self.trace();
        let i2 = self.at(0, 0) * self.at(1, 1)
            + self.at(1, 1) * self.at(2, 2)
            + self.at(2, 2) * self.at(0, 0)
            - (self.at(0, 1) * self.at(1, 0)
                + self.at(1, 2) * self.at(2, 1)
                + self.at(0, 2) * self.at(2, 0));
        let i3 = self.det();

        // Depress the cubic with λ = m + I1/3:  m³ + p m + q = 0
        let m = i1 / 3.0;
        let p = i2 - i1 * i1 / 3.0;
        let q = -2.0 * i1 * i1 * i1 / 27.0 + i1 * i2 / 3.0 - i3;

        // For symmetric matrices p <= 0 and all roots are real:
        //   m_k = 2·sqrt(p_c/3)·cos( (1/3)·arccos( -3√3·q / (2·p_c^{3/2}) ) − 2πk/3 )
        // with p_c = −p. The clamp only guards roundoff on p ≈ 0.
        let p_c = (-p).max(0.0);
        if p_c < crate::EPS {
            return [m, m, m];
        }
        let p32 = p_c * p_c.sqrt();
        let arg = (-3.0f64.sqrt() * 3.0 * q / (2.0 * p32)).clamp(-1.0, 1.0);
        let theta = arg.acos() / 3.0;
        let a = 2.0 * (p_c / 3.0).sqrt();
        let two_pi = 2.0 * core::f64::consts::PI;
        let mut roots = [
            m + a * theta.cos(),
            m + a * (theta - two_pi / 3.0).cos(),
            m + a * (theta - two_pi * 2.0 / 3.0).cos(),
        ];
        roots.sort_by(|x, y| y.partial_cmp(x).unwrap_or(core::cmp::Ordering::Equal));
        roots
    }
}

impl core::ops::Mul<f64> for Mat3 {
    type Output = Mat3;
    #[inline]
    fn mul(self, s: f64) -> Mat3 {
        Mat3::from_rows(self.r0 * s, self.r1 * s, self.r2 * s)
    }
}

impl core::ops::Mul<Mat3> for f64 {
    type Output = Mat3;
    #[inline]
    fn mul(self, m: Mat3) -> Mat3 {
        m * self
    }
}

impl core::ops::Mul<Mat3> for Mat3 {
    type Output = Mat3;
    #[inline]
    fn mul(self, rhs: Mat3) -> Mat3 {
        self.mul_mat(&rhs)
    }
}

impl core::ops::Mul<Vec3> for Mat3 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        self.mul_vec(v)
    }
}

impl core::ops::Add for Mat3 {
    type Output = Mat3;
    #[inline]
    fn add(self, rhs: Mat3) -> Mat3 {
        Mat3::from_rows(self.r0 + rhs.r0, self.r1 + rhs.r1, self.r2 + rhs.r2)
    }
}

impl core::ops::Sub for Mat3 {
    type Output = Mat3;
    #[inline]
    fn sub(self, rhs: Mat3) -> Mat3 {
        Mat3::from_rows(self.r0 - rhs.r0, self.r1 - rhs.r1, self.r2 - rhs.r2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn determinant_and_inverse() {
        let m = Mat3::from_array([4.0, 7.0, 2.0, 3.0, 6.0, 1.0, 2.0, 5.0, 3.0]);
        assert!(close(m.det(), 9.0));
        let inv = m.inverse().expect("invertible");
        let prod = m.mul_mat(&inv);
        for i in 0..3 {
            for j in 0..3 {
                let expect = if i == j { 1.0 } else { 0.0 };
                assert!(close(prod.at(i, j), expect), "({i},{j})={}", prod.at(i, j));
            }
        }
    }

    #[test]
    fn singular_inverse_is_none() {
        let m = Mat3::from_rows(Vec3::new(1.0, 2.0, 3.0), Vec3::new(2.0, 4.0, 6.0), Vec3::X);
        assert!(m.inverse().is_none());
    }

    #[test]
    fn rotation_preserves_norm_and_orthogonality() {
        let m = Mat3::rotation_axis_angle(Vec3::new(1.0, 1.0, 0.0), 0.7);
        let v = Vec3::new(3.0, -2.0, 5.0);
        assert!(close((m * v).norm_squared(), v.norm_squared()));
        let rrt = m.mul_mat(&m.transpose());
        for i in 0..3 {
            for j in 0..3 {
                let expect = if i == j { 1.0 } else { 0.0 };
                assert!((rrt.at(i, j) - expect).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn symmetric_eigenvalues_diagonal() {
        let m = Mat3::diagonal([3.0, -1.0, 7.0]);
        let e = m.symmetric_eigenvalues();
        assert!(close(e[0], 7.0), "{e:?}");
        assert!(close(e[1], 3.0), "{e:?}");
        assert!(close(e[2], -1.0), "{e:?}");
    }

    #[test]
    fn symmetric_eigenvalues_shear_plus_tension() {
        // Eigenvalues of [[2,1,0],[1,2,0],[0,0,5]] are 5, 3, 1.
        let m = Mat3::from_array([2.0, 1.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 5.0]);
        let e = m.symmetric_eigenvalues();
        assert!(close(e[0], 5.0), "{e:?}");
        assert!(close(e[1], 3.0), "{e:?}");
        assert!(close(e[2], 1.0), "{e:?}");
    }

    #[test]
    fn eigenvalue_invariants_hold() {
        let m = Mat3::from_array([4.0, 1.0, 0.5, 1.0, -2.0, 0.3, 0.5, 0.3, 1.0]);
        let e = m.symmetric_eigenvalues();
        let sum: f64 = e.iter().sum();
        let prod: f64 = e.iter().product();
        assert!(close(sum, m.trace()), "trace {sum}");
        assert!(close(prod, m.det()), "det {prod}");
    }
}
