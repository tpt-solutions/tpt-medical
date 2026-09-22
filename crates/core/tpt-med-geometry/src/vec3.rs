//! 3-vector arithmetic.

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

/// A 3-component vector (`f64`), used for positions, directions, and free
/// vectors alike. Lengths follow the caller's unit convention (mm).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// Origin / zero vector.
    pub const ZERO: Vec3 = Vec3::new(0.0, 0.0, 0.0);
    /// Unit +X.
    pub const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
    /// Unit +Y.
    pub const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
    /// Unit +Z.
    pub const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

    /// Constructs a vector from components.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Splat: same value in every component.
    #[inline]
    pub const fn splat(v: f64) -> Self {
        Self::new(v, v, v)
    }

    /// Component-wise access as an array.
    #[inline]
    pub const fn to_array(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    /// Euclidean norm.
    #[inline]
    pub fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    /// Squared norm (avoids the sqrt in comparisons).
    #[inline]
    pub fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product (`self × o`), right-handed.
    #[inline]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Unit-length copy; the zero vector is returned unchanged.
    #[inline]
    pub fn normalize(self) -> Self {
        let n = self.norm();
        if n > crate::EPS {
            self / n
        } else {
            self
        }
    }

    /// Linear interpolation `self + t * (other - self)`.
    #[inline]
    pub fn lerp(self, other: Self, t: f64) -> Self {
        self + (other - self) * t
    }

    /// Distance to another point.
    #[inline]
    pub fn distance(self, other: Self) -> f64 {
        (other - self).norm()
    }

    /// Component-wise minimum.
    #[inline]
    pub fn min_component(self) -> f64 {
        self.x.min(self.y).min(self.z)
    }

    /// Component-wise maximum.
    #[inline]
    pub fn max_component(self) -> f64 {
        self.x.max(self.y).max(self.z)
    }

    /// Component-wise product with another vector.
    #[inline]
    pub fn component_mul(self, o: Self) -> Self {
        Self::new(self.x * o.x, self.y * o.y, self.z * o.z)
    }

    /// Component-wise division by another vector.
    ///
    /// Panics only through IEEE on zero components (inf/NaN), matching the
    /// rest of the numeric stack.
    #[inline]
    pub fn component_div(self, o: Self) -> Self {
        Self::new(self.x / o.x, self.y / o.y, self.z / o.z)
    }
}

impl Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn div(self, s: f64) -> Self {
        Self::new(self.x / s, self.y / s, self.z / s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    #[test]
    fn cross_is_anticommutative_and_right_handed() {
        let c = Vec3::X.cross(Vec3::Y);
        assert_eq!(c, Vec3::Z);
        assert_eq!(Vec3::Y.cross(Vec3::X), -Vec3::Z);
    }

    #[test]
    fn norm_and_normalize() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert!(close(v.norm(), 5.0));
        let u = v.normalize();
        assert!(close(u.norm(), 1.0));
        assert_eq!(Vec3::ZERO.normalize(), Vec3::ZERO);
    }

    #[test]
    fn lerp_endpoints() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(10.0, 20.0, 30.0);
        assert_eq!(a.lerp(b, 0.0), a);
        assert_eq!(a.lerp(b, 1.0), b);
        assert_eq!(a.lerp(b, 0.5), Vec3::new(5.0, 10.0, 15.0));
    }
}
