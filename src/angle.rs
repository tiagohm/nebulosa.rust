use std::{
    f64::consts::PI,
    ops::{Div, Mul},
    process::Output,
};

// Radians to degrees.
pub const RAD2DEG: f64 = 180.0 / PI;

// Degrees to radians.
pub const DEG2RAD: f64 = PI / 180.0;

// Arcminutes to radians.
pub const AMIN2RAD: f64 = PI / 180.0 / 60.0;

// Arcsecconds to radians.
pub const ASEC2RAD: f64 = PI / 180.0 / 3600.0;

// Milliarcsecconds to radians.
pub const MILLIASEC2RAD: f64 = PI / 180.0 / 3600000.0;

// Angular velocity in radians/s.
pub const ANGULAR_VELOCITY: f64 = 7.2921150E-5;

// Arcseconds in a full circle.
pub const TURNAS: f64 = 1296000.0;

pub trait Angle {
    fn radians(&self) -> Radians;

    fn degress(&self) -> Degress;

    fn arcmin(&self) -> ArcMin;

    fn arcsec(&self) -> ArcSec;

    fn mas(&self) -> Mas;
}

#[derive(Default, Clone, Copy)]
pub struct Radians(pub f64);
#[derive(Default, Clone, Copy)]
pub struct Degress(pub f64);
#[derive(Default, Clone, Copy)]
pub struct ArcMin(pub f64);
#[derive(Default, Clone, Copy)]
pub struct ArcSec(pub f64);
#[derive(Default, Clone, Copy)]
pub struct Mas(pub f64);

impl Radians {
    #[inline(always)]
    fn new(value: f64) -> Radians {
        Radians(value)
    }

    #[inline(always)]
    fn degrees(value: f64) -> Radians {
        Radians(value * DEG2RAD)
    }

    #[inline(always)]
    fn arcmin(value: f64) -> Radians {
        Radians(value * AMIN2RAD)
    }

    #[inline(always)]
    fn arcsec(value: f64) -> Radians {
        Radians(value * ASEC2RAD)
    }

    #[inline(always)]
    fn mas(value: f64) -> Radians {
        Radians(value * MILLIASEC2RAD)
    }
}

impl Angle for Radians {
    #[inline(always)]
    fn radians(&self) -> Radians {
        *self
    }

    #[inline(always)]
    fn degress(&self) -> Degress {
        Degress(self.0 * RAD2DEG)
    }

    #[inline(always)]
    fn arcmin(&self) -> ArcMin {
        ArcMin(self.0 * (60.0 * RAD2DEG))
    }

    #[inline(always)]
    fn arcsec(&self) -> ArcSec {
        ArcSec(self.0 * (3600.0 * RAD2DEG))
    }

    #[inline(always)]
    fn mas(&self) -> Mas {
        Mas(self.0 * (3600000.0 * RAD2DEG))
    }
}

impl Mul<f64> for Radians {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        Radians(self.0 * rhs)
    }
}

impl Div<f64> for Radians {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        Radians(self.0 / rhs)
    }
}

impl ArcSec {
    #[inline(always)]
    fn new(value: f64) -> ArcSec {
        ArcSec(value)
    }

    #[inline(always)]
    fn degrees(value: f64) -> ArcSec {
        ArcSec(value * 3600.0)
    }

    #[inline(always)]
    fn radians(value: f64) -> ArcSec {
        ArcSec(value / ASEC2RAD)
    }

    #[inline(always)]
    fn arcmin(value: f64) -> ArcSec {
        ArcSec(value * 60.0)
    }

    #[inline(always)]
    fn mas(value: f64) -> ArcSec {
        ArcSec(value / 1000.0)
    }
}

impl Angle for ArcSec {
    #[inline(always)]
    fn radians(&self) -> Radians {
        Radians::arcsec(self.0)
    }

    #[inline(always)]
    fn degress(&self) -> Degress {
        todo!()
    }

    #[inline(always)]
    fn arcmin(&self) -> ArcMin {
        todo!()
    }

    #[inline(always)]
    fn arcsec(&self) -> ArcSec {
        *self
    }

    #[inline(always)]
    fn mas(&self) -> Mas {
        todo!()
    }
}

impl Mul<f64> for ArcSec {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self::Output {
        ArcSec(self.0 * rhs)
    }
}

impl Div<f64> for ArcSec {
    type Output = Self;

    fn div(self, rhs: f64) -> Self::Output {
        ArcSec(self.0 / rhs)
    }
}

#[cfg(test)]
mod test {
    use assertor::{assert_that, EqualityAssertion, FloatAssertion};

    use crate::angle::Radians;

    use super::Angle;

    #[test]
    fn radians() {
        assert_that!(Radians::new(1.4).0).is_equal_to(1.4);
        assert_that!(Radians::degrees(45.0).0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.7853981633974);
        assert_that!(Radians::arcmin(90.0).0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.026179938779914);
        assert_that!(Radians::arcsec(1800.0).0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.0087266462599716);
        assert_that!(Radians::mas(4500000.0).0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.02181661564992911);
        assert_that!((Radians::new(1.4) * 2.0).0).is_equal_to(2.8);
        assert_that!((Radians::new(1.4) / 1.4).0).is_equal_to(1.0);

        assert_that!(Radians(1.4).0).is_equal_to(1.4);
        assert_that!(Radians(1.4).degress().0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(80.2140913183152);
        assert_that!(Radians(1.4).arcmin().0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(4812.845479098914953);
        assert_that!(Radians(1.4).arcsec().0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(288770.728745934897219);
        assert_that!(Radians(1.4).mas().0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(288770728.745934897219);
    }
}
