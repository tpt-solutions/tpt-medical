//! Type-safe unit system for biomedical simulation.
//!
//! Biomechanics mixes SI mechanics, scanner conventions, and physiology
//! conventions (mmHg). This crate provides cheap newtype wrappers over `f64`
//! for the quantities that flow through the stack, with explicit conversions
//! and no implicit cross-unit arithmetic. Base conventions:
//!
//! | Quantity      | Canonical unit |
//! |---------------|----------------|
//! | Length        | millimetre     |
//! | Pressure      | megapascal     |
//! | Force         | newton         |
//! | Density       | g/cm³          |
//! | Time          | second         |
//! | Angle         | radian         |
//! | Viscosity     | Pa·s           |
//!
//! # Examples
//!
//! ```
//! use tpt_med_units::{Length, Pressure};
//!
//! let lesion = Length::from_mm(30.0);
//! let systolic = Pressure::from_mmhg(120.0);
//! assert!((systolic.to_pascal() - 15_998.7).abs() < 1.0);
//! assert_eq!(lesion.to_cm(), 3.0);
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use core::fmt;
use core::ops::{Add, Div, Mul, Sub};

macro_rules! quantity {
    ($(#[$meta:meta])* $name:ident, $unit:expr, $from:ident, $to:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, Default, PartialEq, PartialOrd)]
        #[repr(transparent)]
        pub struct $name(pub f64);

        impl $name {
            /// Zero-valued quantity.
            pub const ZERO: Self = Self(0.0);

            /// Creates a quantity from the canonical unit value.
            #[inline]
            pub const fn new(value: f64) -> Self {
                Self(value)
            }

            /// Creates a quantity from the canonical unit value
            /// (alias of [`Self::new`], reads naturally at call sites).
            #[inline]
            pub const fn $from(value: f64) -> Self {
                Self(value)
            }

            /// Raw value in the canonical unit (alias of [`Self::value`]).
            #[inline]
            pub const fn $to(self) -> f64 {
                self.0
            }

            /// Raw value in the canonical unit.
            #[inline]
            pub const fn value(self) -> f64 {
                self.0
            }

            /// Canonical unit symbol (for display only).
            pub const UNIT: &'static str = $unit;

            /// Absolute value.
            #[inline]
            pub fn abs(self) -> Self {
                Self(self.0.abs())
            }

            /// Largest of two quantities.
            #[inline]
            pub fn max(self, other: Self) -> Self {
                Self(self.0.max(other.0))
            }

            /// Smallest of two quantities.
            #[inline]
            pub fn min(self, other: Self) -> Self {
                Self(self.0.min(other.0))
            }

            /// True if the value is finite (not NaN or infinite).
            #[inline]
            pub fn is_finite(self) -> bool {
                self.0.is_finite()
            }
        }

        impl Add for $name {
            type Output = Self;
            #[inline]
            fn add(self, rhs: Self) -> Self {
                Self(self.0 + rhs.0)
            }
        }

        impl Sub for $name {
            type Output = Self;
            #[inline]
            fn sub(self, rhs: Self) -> Self {
                Self(self.0 - rhs.0)
            }
        }

        impl Mul<f64> for $name {
            type Output = Self;
            #[inline]
            fn mul(self, rhs: f64) -> Self {
                Self(self.0 * rhs)
            }
        }

        impl Div<f64> for $name {
            type Output = Self;
            #[inline]
            fn div(self, rhs: f64) -> Self {
                Self(self.0 / rhs)
            }
        }

        impl Div for $name {
            type Output = f64;
            #[inline]
            fn div(self, rhs: Self) -> f64 {
                self.0 / rhs.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{:.6} {}", self.0, $unit)
            }
        }
    };
}

quantity!(
    /// Length in millimetres.
    Length,
    "mm",
    from_mm,
    to_mm
);
quantity!(
    /// Pressure in megapascals.
    Pressure,
    "MPa",
    from_mpa,
    to_mpa
);
quantity!(
    /// Force in newtons.
    Force,
    "N",
    from_n,
    to_n
);
quantity!(
    /// Mass density in g/cm³.
    Density,
    "g/cm3",
    from_gcm3,
    to_gcm3
);
quantity!(
    /// Time in seconds.
    Time,
    "s",
    from_s,
    to_s
);
quantity!(
    /// Angle in radians.
    Angle,
    "rad",
    from_rad,
    to_rad
);
quantity!(
    /// Dynamic viscosity in Pa·s.
    Viscosity,
    "Pa.s",
    from_pas,
    to_pas
);
quantity!(
    /// Velocity in mm/s.
    Velocity,
    "mm/s",
    from_mms,
    to_mms
);
quantity!(
    /// Young's modulus in megapascals.
    Modulus,
    "MPa",
    from_mpa,
    to_mpa
);
quantity!(
    /// Volumetric flow rate in mm³/s.
    FlowRate,
    "mm3/s",
    from_mm3s,
    to_mm3s
);

impl Length {
    /// Creates a length from centimetres.
    #[inline]
    pub const fn from_cm(cm: f64) -> Self {
        Self(cm * 10.0)
    }

    /// Converts to centimetres.
    #[inline]
    pub const fn to_cm(self) -> f64 {
        self.0 / 10.0
    }

    /// Creates a length from metres.
    #[inline]
    pub const fn from_m(m: f64) -> Self {
        Self(m * 1000.0)
    }

    /// Converts to metres.
    #[inline]
    pub const fn to_m(self) -> f64 {
        self.0 / 1000.0
    }
}

impl Pressure {
    /// Creates a pressure from pascals.
    #[inline]
    pub const fn from_pa(pa: f64) -> Self {
        Self(pa / 1.0e6)
    }

    /// Converts to pascals.
    #[inline]
    pub const fn to_pascal(self) -> f64 {
        self.0 * 1.0e6
    }

    /// Creates a pressure from mmHg (1 mmHg = 133.322387415 Pa, the
    /// conventional millimetre of mercury).
    #[inline]
    pub fn from_mmhg(mmhg: f64) -> Self {
        Self(mmhg * 133.322_387_415 / 1.0e6)
    }

    /// Converts to mmHg.
    #[inline]
    pub fn to_mmhg(self) -> f64 {
        self.0 * 1.0e6 / 133.322_387_415
    }

    /// Creates a pressure from kilopascals.
    #[inline]
    pub const fn from_kpa(kpa: f64) -> Self {
        Self(kpa / 1000.0)
    }

    /// Converts to kilopascals.
    #[inline]
    pub const fn to_kpa(self) -> f64 {
        self.0 * 1000.0
    }
}

impl Force {
    /// Creates a force from kilonewtons (joint-reaction load tables).
    #[inline]
    pub const fn from_kn(kn: f64) -> Self {
        Self(kn * 1000.0)
    }

    /// Converts to kilonewtons.
    #[inline]
    pub const fn to_kn(self) -> f64 {
        self.0 / 1000.0
    }

    /// Force = multiplier × body weight (gait-load convention).
    pub fn body_weights(bw_multiple: f64, body_weight: Self) -> Self {
        Self(bw_multiple * body_weight.0)
    }
}

impl Time {
    /// Creates a time from milliseconds (ECG/gait timing conventions).
    #[inline]
    pub const fn from_ms(ms: f64) -> Self {
        Self(ms / 1000.0)
    }

    /// Converts to milliseconds.
    #[inline]
    pub const fn to_ms(self) -> f64 {
        self.0 * 1000.0
    }

    /// Creates a time from minutes (heart-rate period conventions).
    #[inline]
    pub fn from_min(min: f64) -> Self {
        Self(min * 60.0)
    }
}

impl Viscosity {
    /// Creates a viscosity from centipoise (blood rheology convention;
    /// 1 cP = 1 mPa·s).
    #[inline]
    pub const fn from_cp(cp: f64) -> Self {
        Self(cp / 1000.0)
    }

    /// Converts to centipoise.
    #[inline]
    pub const fn to_cp(self) -> f64 {
        self.0 * 1000.0
    }
}

impl FlowRate {
    /// Creates a flow rate from mL/min (cardiac catheterisation convention).
    #[inline]
    pub fn from_ml_per_min(ml_min: f64) -> Self {
        // 1 mL = 1000 mm³, 1 min = 60 s
        Self(ml_min * 1000.0 / 60.0)
    }

    /// Converts to mL/min.
    #[inline]
    pub fn to_ml_per_min(self) -> f64 {
        self.0 * 60.0 / 1000.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol * b.abs().max(1.0)
    }

    #[test]
    fn pressure_conversions_roundtrip() {
        let p = Pressure::from_mmhg(120.0);
        assert!(close(p.to_pascal(), 15_998.686_5, 1e-4));
        assert!(close(p.to_mmhg(), 120.0, 1e-12));
        assert!(close(Pressure::from_pa(101_325.0).to_kpa(), 101.325, 1e-9));
    }

    #[test]
    fn length_conversions() {
        assert_eq!(Length::from_cm(2.5).to_mm(), 25.0);
        assert_eq!(Length::from_m(1.0).to_cm(), 100.0);
    }

    #[test]
    fn force_body_weight() {
        // Typical stance peak: ~3×BW for a 75 kg patient (BW ≈ 735.75 N).
        let bw = Force::from_n(75.0 * 9.81);
        let joint = Force::body_weights(3.0, bw);
        assert!(close(joint.to_kn(), 2.207_25, 1e-4));
    }

    #[test]
    fn time_and_viscosity() {
        assert_eq!(Time::from_ms(800.0).to_s(), 0.8);
        assert_eq!(Time::from_min(1.0).to_s(), 60.0);
        // Whole blood ~3.5 cP at high shear
        assert!(close(Viscosity::from_cp(3.5).to_pas(), 0.0035, 1e-12));
    }

    #[test]
    fn flow_rate_ml_per_min() {
        // 1000 mL/min = 1e6 mm³ / 60 s ≈ 16666.67 mm³/s
        let q = FlowRate::from_ml_per_min(1000.0);
        assert!(close(q.to_mm3s(), 16_666.666_7, 1e-4));
        assert!(close(q.to_ml_per_min(), 1000.0, 1e-9));
    }

    #[test]
    fn arithmetic_operators() {
        let a = Length::from_mm(3.0);
        let b = Length::from_mm(2.0);
        assert_eq!((a + b).to_mm(), 5.0);
        assert_eq!((a - b).to_mm(), 1.0);
        assert_eq!((a * 2.0).to_mm(), 6.0);
        assert_eq!(a / b, 1.5);
        assert_eq!(a.max(b), a);
        assert_eq!(a.min(b), b);
    }

    #[test]
    fn display_includes_unit() {
        let s = format!("{}", Pressure::from_mpa(1.5));
        assert!(s.contains("MPa"), "got {s}");
    }
}
