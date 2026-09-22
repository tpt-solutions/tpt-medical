//! Axis-aligned bounding box.

use crate::vec3::Vec3;

/// Axis-aligned bounding box, inclusive of its bounds. Used for spatial
/// pruning, mesh extents, and viewport fitting.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Aabb {
    /// Minimum corner.
    pub min: Vec3,
    /// Maximum corner.
    pub max: Vec3,
}

impl Aabb {
    /// Degenerate box around a point.
    pub fn point(p: Vec3) -> Self {
        Self { min: p, max: p }
    }

    /// Smallest box containing all points. Empty input yields a degenerate
    /// box at the origin.
    pub fn from_points<'a>(points: impl IntoIterator<Item = &'a Vec3>) -> Self {
        let mut iter = points.into_iter();
        match iter.next() {
            None => Self::point(Vec3::ZERO),
            Some(first) => {
                let mut b = Self::point(*first);
                for p in iter {
                    b.expand(*p);
                }
                b
            }
        }
    }

    /// Grows the box to contain a point.
    pub fn expand(&mut self, p: Vec3) {
        self.min = Vec3::new(
            self.min.x.min(p.x),
            self.min.y.min(p.y),
            self.min.z.min(p.z),
        );
        self.max = Vec3::new(
            self.max.x.max(p.x),
            self.max.y.max(p.y),
            self.max.z.max(p.z),
        );
    }

    /// Union with another box.
    pub fn union(&self, other: &Aabb) -> Aabb {
        let mut b = *self;
        b.expand(other.min);
        b.expand(other.max);
        b
    }

    /// Center of the box.
    pub fn center(&self) -> Vec3 {
        (self.min + self.max) * 0.5
    }

    /// Extents (full side lengths).
    pub fn extents(&self) -> Vec3 {
        self.max - self.min
    }

    /// True if the point is inside (inclusive).
    pub fn contains(&self, p: Vec3) -> bool {
        p.x >= self.min.x
            && p.y >= self.min.y
            && p.z >= self.min.z
            && p.x <= self.max.x
            && p.y <= self.max.y
            && p.z <= self.max.z
    }

    /// True if the boxes overlap (inclusive of touching faces).
    pub fn intersects(&self, other: &Aabb) -> bool {
        !(self.max.x < other.min.x
            || self.min.x > other.max.x
            || self.max.y < other.min.y
            || self.min.y > other.max.y
            || self.max.z < other.min.z
            || self.min.z > other.max.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_points_and_expand() {
        let pts = [
            Vec3::new(1.0, -2.0, 3.0),
            Vec3::new(-4.0, 5.0, 0.0),
            Vec3::new(0.0, 0.0, 7.0),
        ];
        let b = Aabb::from_points(pts.iter());
        assert_eq!(b.min, Vec3::new(-4.0, -2.0, 0.0));
        assert_eq!(b.max, Vec3::new(1.0, 5.0, 7.0));
        assert_eq!(b.center(), Vec3::new(-1.5, 1.5, 3.5));
        assert_eq!(b.extents(), Vec3::new(5.0, 7.0, 7.0));
    }

    #[test]
    fn contains_and_intersects() {
        let a = Aabb {
            min: Vec3::ZERO,
            max: Vec3::splat(10.0),
        };
        assert!(a.contains(Vec3::new(5.0, 5.0, 5.0)));
        assert!(a.contains(a.min));
        assert!(!a.contains(Vec3::new(10.1, 0.0, 0.0)));
        let b = Aabb {
            min: Vec3::new(9.0, 9.0, 9.0),
            max: Vec3::splat(20.0),
        };
        assert!(a.intersects(&b));
        let c = Aabb {
            min: Vec3::new(11.0, 11.0, 11.0),
            max: Vec3::splat(20.0),
        };
        assert!(!a.intersects(&c));
    }
}
