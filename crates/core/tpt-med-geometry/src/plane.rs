//! Infinite plane primitive.

use crate::vec3::Vec3;

/// A plane defined by a unit normal `n` and offset `d`, such that points `p`
/// with `n·p + d = 0` lie on the plane. `d` is the signed distance from the
/// origin along the normal. Used for osteotomy cuts and clipping.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    /// Unit normal.
    pub normal: Vec3,
    /// Signed offset from origin along the normal.
    pub offset: f64,
}

impl Plane {
    /// Builds a plane from a (non-zero) normal and one point on the plane.
    /// The normal is normalized; returns `None` for a zero normal.
    pub fn from_point_normal(point: Vec3, normal: Vec3) -> Option<Self> {
        let n = normal.normalize();
        if n.norm_squared() < crate::EPS {
            return None;
        }
        Some(Self {
            normal: n,
            offset: -n.dot(point),
        })
    }

    /// Builds a plane through three points (counter-clockwise when viewed
    /// from the side the normal points to). Returns `None` for degenerate
    /// (collinear) triples.
    pub fn from_three_points(a: Vec3, b: Vec3, c: Vec3) -> Option<Self> {
        let normal = (b - a).cross(c - a);
        Self::from_point_normal(a, normal)
    }

    /// Signed distance from a point: positive on the normal side.
    #[inline]
    pub fn signed_distance(&self, p: Vec3) -> f64 {
        self.normal.dot(p) + self.offset
    }

    /// Ray–plane intersection: `origin + t · dir` where the ray hits the
    /// plane. Returns `None` when the ray is parallel or points away.
    pub fn ray_intersection(&self, origin: Vec3, dir: Vec3) -> Option<(Vec3, f64)> {
        let denom = self.normal.dot(dir);
        if denom.abs() < crate::EPS {
            return None;
        }
        let t = -self.signed_distance(origin) / denom;
        if t < 0.0 {
            return None;
        }
        Some((origin + dir * t, t))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signed_distance_signs() {
        // z = 0 plane with +Z normal
        let plane = Plane::from_point_normal(Vec3::ZERO, Vec3::Z).unwrap();
        assert_eq!(plane.signed_distance(Vec3::new(5.0, -3.0, 2.0)), 2.0);
        assert_eq!(plane.signed_distance(Vec3::new(5.0, -3.0, -2.0)), -2.0);
        assert_eq!(plane.signed_distance(Vec3::new(5.0, -3.0, 0.0)), 0.0);
    }

    #[test]
    fn offset_plane() {
        // z = 4 plane with +Z normal: offset = -4
        let plane = Plane::from_point_normal(Vec3::new(0.0, 0.0, 4.0), Vec3::Z).unwrap();
        assert!((plane.offset + 4.0).abs() < 1e-12);
        assert_eq!(plane.signed_distance(Vec3::new(0.0, 0.0, 6.0)), 2.0);
    }

    #[test]
    fn three_points_degenerate_is_none() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let c = Vec3::new(2.0, 0.0, 0.0);
        assert!(Plane::from_three_points(a, b, c).is_none());
        assert!(Plane::from_three_points(a, b, Vec3::Y).is_some());
    }

    #[test]
    fn ray_hits_plane() {
        let plane = Plane::from_point_normal(Vec3::new(0.0, 0.0, 10.0), -Vec3::Z).unwrap();
        let (hit, t) = plane
            .ray_intersection(Vec3::new(1.0, 1.0, 0.0), Vec3::Z)
            .expect("hits");
        assert_eq!(hit, Vec3::new(1.0, 1.0, 10.0));
        assert!((t - 10.0).abs() < 1e-12);
        // Parallel ray misses
        assert!(plane.ray_intersection(Vec3::ZERO, Vec3::X).is_none());
    }
}
