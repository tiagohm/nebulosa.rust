use std::cmp::{max, min};

use crate::{
    angle::{ArcSec, Radians},
    time::{DAYSEC, DAYSPERJC, J2000, MJD0},
};

const DBL_EPSILON: f64 = 2.220446049250313E-16;

pub type LeapSecondChange = (i32, i32, f64);
pub type LeapSecondDrift = (f64, f64);

pub fn round_to_nearest_whole_number(a: f64) -> f64 {
    if (f64::abs(a) < 0.5) {
        0.0
    } else if a < 0.0 {
        f64::ceil(a - 0.5)
    } else {
        f64::floor(a + 0.5)
    }
}

/// International Atomic Time, TAI, to Universal Time, UT1.
#[inline(always)]
pub fn era_tai_ut1(tai1: f64, tai2: f64, ut1_minus_tai: f64) -> (f64, f64) {
    (tai1, tai2 + ut1_minus_tai / DAYSEC)
}

/// Universal Time, UT1, to International Atomic Time, TAI.
#[inline(always)]
pub fn era_ut1_tai(ut11: f64, ut12: f64, ut1_minus_tai: f64) -> (f64, f64) {
    (ut11, ut12 - ut1_minus_tai / DAYSEC)
}

/// International Atomic Time, TAI, to Coordinated Universal Time, UTC.
pub fn era_tai_utc(tai1: f64, tai2: f64) -> (f64, f64) {
    let mut u2 = tai2;

    // Iterate(though in most cases just once is enough).
    for i in 0..3 {
        let (g1, g2) = era_utc_tai(tai1, u2);

        // Adjust guessed UTC.
        u2 += tai1 - g1;
        u2 += tai2 - g2;
    }

    (tai1, u2)
}

pub fn era_utc_tai(utc1: f64, utc2: f64) -> (f64, f64) {
    let u1 = f64::max(utc1, utc2);
    let u2 = f64::min(utc1, utc2);

    // Get TAI-UTC at 0h today.
    let cal = era_jd_to_cal(u1, u2);
    let dat0 = era_dat(cal.0, cal.1, cal.2, 0.0);

    // Get TAI-UTC at 12h today (to detect drift).
    let dat12 = era_dat(cal.0, cal.1, cal.2, 0.5);

    // Get TAI-UTC at 0h tomorrow (to detect jumps).
    let calt = era_jd_to_cal(u1 + 1.5, u2 - cal.3);
    let dat24 = era_dat(cal.0, cal.1, cal.2, 0.0);

    // Separate TAI-UTC change into per-day (DLOD) and any jump (DLEAP).
    let dlod = 2.0 * (dat12 - dat0);
    let dleap = dat24 - (dat0 + dlod);

    // Remove any scaling applied to spread leap into preceding day.
    let mut fd = cal.3 * (DAYSEC + dleap) / DAYSEC;

    // Scale from (pre-1972) UTC seconds to SI seconds.
    fd *= (DAYSEC + dlod) / DAYSEC;

    // Today's calendar date to 2-part JD.
    let z = era_cal_to_jd(cal.0, cal.1, cal.2);

    // Assemble the TAI result, preserving the UTC split and order.
    let a2 = (MJD0 - u1) + z + (fd + dat0 / DAYSEC);

    (u1, a2)
}

pub fn era_utc_ut1(utc1: f64, utc2: f64, dut1: f64) -> (f64, f64) {
    let cal = era_jd_to_cal(utc1, utc2);
    let dat = era_dat(cal.0, cal.1, cal.2, cal.3);

    // Form UT1-TAI
    let dta = dut1 - dat;

    let (tai1, tai2) = era_utc_tai(utc1, utc2);
    era_tai_ut1(tai1, tai2, dta)
}

pub fn era_ut1_utc(ut11: f64, ut12: f64, dut1: f64) -> (f64, f64) {
    let u1 = f64::max(ut11, ut12);
    let mut u2 = f64::min(ut11, ut12);

    let mut duts = dut1;

    // See if the UT1 can possibly be in a leap-second day.
    let mut d1 = u1;
    let mut dats1 = 0.0;

    for i in -1..4 {
        let mut d2 = u2 + i as f64;
        let cal = era_jd_to_cal(d1, d2);
        let dats2 = era_dat(cal.0, cal.1, cal.2, 0.0);

        if i == -1 {
            dats1 = dats2
        }

        let ddats = dats2 - dats1;

        if (f64::abs(ddats) >= 0.5) {
            // Yes, leap second nearby: ensure UT1-UTC is "before" value.
            if ddats * duts >= 0.0 {
                duts -= ddats;
            }

            // UT1 for the start of the UTC day that ends in a leap.
            d1 = MJD0;
            d2 = era_cal_to_jd(cal.0, cal.1, cal.2);

            let us1 = d1;
            let us2 = d2 - 1.0 + duts / DAYSEC;

            // Is the UT1 after this point?
            let du = u1 - us1 + (u2 - us2);

            if du > 0.0 {
                // Yes: fraction of the current UTC day that has elapsed.
                let fd = du * DAYSEC / (DAYSEC + ddats);

                // Ramp UT1-UTC to bring about ERFA's JD(UTC) convention.
                duts += ddats * if (fd <= 1.0) { fd } else { 1.0 };
            }

            break;
        }

        dats1 = dats2
    }

    // Subtract the (possibly adjusted) UT1-UTC from UT1 to give UTC.
    u2 -= duts / DAYSEC;

    (u1, u2)
}

pub fn era_jd_to_cal(dj1: f64, dj2: f64) -> (i32, i32, i32, f64) {
    // Separate day and fraction (where -0.5 <= fraction < 0.5).
    let d = round_to_nearest_whole_number(dj1);

    let f1 = dj1 - d;
    let mut jd = d as i64;
    let d = round_to_nearest_whole_number(dj2);
    let f2 = dj2 - d;
    jd += d as i64;

    // Compute f1+f2+0.5 using compensated summation (Klein 2006).
    let mut s = 0.5;
    let mut cs = 0.0;
    let v = [f1, f2];

    for i in 0..2 {
        let x = v[i];
        let t = s + x;

        cs += if f64::abs(s) >= f64::abs(x) {
            (s - t) + x
        } else {
            (x - t) + s
        };

        s = t;

        if s >= 1.0 {
            jd += 1;
            s -= 1.0;
        }
    }

    let mut f = s + cs;
    cs = f - s;

    // Deal with negative f.
    if f < 0.0 {
        // Compensated summation: assume that |s| <= 1.0.
        f = s + 1.0;
        cs += (1.0 - f) + s;
        s = f;
        f = s + cs;
        cs = f - s;
        jd -= 1;
    }

    // Deal with f that is 1.0 or more (when rounded to double).
    if (f - 1.0) >= -DBL_EPSILON / 4.0 {
        // Compensated summation: assume that |s| <= 1.0.
        let t = s - 1.0;
        cs += (s - t) - 1.0;
        s = t;
        f = s + cs;

        if -DBL_EPSILON / 2.0 < f {
            jd += 1;
            f = f64::max(f, 0.0);
        }
    }

    // Express day in Gregorian calendar.
    let mut l = jd + 68569;
    let n = (4 * l) / 146097;
    l -= (146097 * n + 3) / 4;
    let i = (4000 * (l + 1)) / 1461001;
    l -= (1461 * i) / 4 - 31;
    let k = (80 * l) / 2447;
    let id = (l - (2447 * k) / 80) as i32;
    l = k / 11;
    let im = (k + 2 - 12 * l) as i32;
    let iy = (100 * (n - 49) + i + l) as i32;

    (iy, im, id, f)
}

pub fn era_cal_to_jd(iy: i32, im: i32, id: i32) -> f64 {
    let my = (im - 14) / 12;
    let iypmy = iy + my;
    ((1461 * (iypmy + 4800)) / 4 + (367 * (im - 2 - 12 * my)) / 12
        - (3 * ((iypmy + 4900) / 100)) / 4
        + id
        - 2432076) as f64
}

/// For a given UTC date, calculate Delta(AT) = TAI-UTC.
pub fn era_dat(iy: i32, im: i32, id: i32, fd: f64) -> f64 {
    let djm = era_cal_to_jd(iy, im, id);

    // Combine year and month to form a date-ordered integer...
    let m = 12 * iy + im as i32;

    // ...and use it to find the preceding table entry.
    let i = match LEAP_SECOND_CHANGES
        .iter()
        .rposition(|&x| m >= 12 * x.0 + x.1)
    {
        Some(i) => i,
        None => return f64::NAN,
    };

    // Get the Delta(AT).
    let mut da = LEAP_SECOND_CHANGES[i].2;

    // If pre-1972, adjust for drift.
    if LEAP_SECOND_CHANGES[i].0 < 1972 {
        da += (djm + fd - LEAP_SECOND_DRIFT[i].0) * LEAP_SECOND_DRIFT[i].1
    }

    da
}

/// The TIO locator s', positioning the Terrestrial Intermediate Origin
/// on the equator of the Celestial Intermediate Pole.
pub fn era_sp00(tt1: f64, tt2: f64) -> ArcSec {
    let t = (tt1 - J2000 + tt2) / DAYSPERJC;
    let sp = -47e-6 * t;
    ArcSec(sp)
}

const LEAP_SECOND_CHANGES: [LeapSecondChange; 42] = [
    (1960, 1, 1.4178180),
    (1961, 1, 1.4228180),
    (1961, 8, 1.3728180),
    (1962, 1, 1.8458580),
    (1963, 11, 1.9458580),
    (1964, 1, 3.2401300),
    (1964, 4, 3.3401300),
    (1964, 9, 3.4401300),
    (1965, 1, 3.5401300),
    (1965, 3, 3.6401300),
    (1965, 7, 3.7401300),
    (1965, 9, 3.8401300),
    (1966, 1, 4.3131700),
    (1968, 2, 4.2131700),
    (1972, 1, 10.0),
    (1972, 7, 11.0),
    (1973, 1, 12.0),
    (1974, 1, 13.0),
    (1975, 1, 14.0),
    (1976, 1, 15.0),
    (1977, 1, 16.0),
    (1978, 1, 17.0),
    (1979, 1, 18.0),
    (1980, 1, 19.0),
    (1981, 7, 20.0),
    (1982, 7, 21.0),
    (1983, 7, 22.0),
    (1985, 7, 23.0),
    (1988, 1, 24.0),
    (1990, 1, 25.0),
    (1991, 1, 26.0),
    (1992, 7, 27.0),
    (1993, 7, 28.0),
    (1994, 7, 29.0),
    (1996, 1, 30.0),
    (1997, 7, 31.0),
    (1999, 1, 32.0),
    (2006, 1, 33.0),
    (2009, 1, 34.0),
    (2012, 7, 35.0),
    (2015, 7, 36.0),
    (2017, 1, 37.0),
];

const LEAP_SECOND_DRIFT: [LeapSecondDrift; 14] = [
    (37300.0, 0.0012960),
    (37300.0, 0.0012960),
    (37300.0, 0.0012960),
    (37665.0, 0.0011232),
    (37665.0, 0.0011232),
    (38761.0, 0.0012960),
    (38761.0, 0.0012960),
    (38761.0, 0.0012960),
    (38761.0, 0.0012960),
    (38761.0, 0.0012960),
    (38761.0, 0.0012960),
    (38761.0, 0.0012960),
    (39126.0, 0.0025920),
    (39126.0, 0.0025920),
];

#[cfg(test)]
mod test {
    use crate::{angle::Angle, erfa::era_sp00};

    use super::{
        era_cal_to_jd, era_dat, era_jd_to_cal, era_tai_ut1, era_tai_utc, era_ut1_utc, era_utc_tai,
        era_utc_ut1,
    };
    use assertor::{assert_that, EqualityAssertion, FloatAssertion};

    #[test]
    fn tai_ut1() {
        let (u1, u2) = era_tai_ut1(2453750.5, 0.892482639, -32.6659);
        assert_that!(u1).is_equal_to(2453750.5);
        assert_that!(u2)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.8921045614537037037);
    }

    #[test]
    fn tai_utc() {
        let (u1, u2) = era_tai_utc(2453750.5, 0.892482639);
        assert_that!(u1).is_equal_to(2453750.5);
        assert_that!(u2)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.8921006945555555556);
    }

    #[test]
    fn utc_tai() {
        let (u1, u2) = era_utc_tai(2453750.5, 0.892100694);
        assert_that!(u1).is_equal_to(2453750.5);
        assert_that!(u2)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.8924826384444444444);
    }

    #[test]
    fn utc_ut1() {
        let (u1, u2) = era_utc_ut1(2453750.5, 0.892100694, 0.3341);
        assert_that!(u1).is_equal_to(2453750.5);
        assert_that!(u2)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.8921045608981481481);
    }

    #[test]
    fn ut1_utc() {
        let (u1, u2) = era_ut1_utc(2453750.5, 0.892104561, 0.3341);
        assert_that!(u1).is_equal_to(2453750.5);
        assert_that!(u2)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(0.8921006941018518519);
    }

    #[test]
    fn dat() {
        assert_that!(era_dat(2003, 6, 1, 0.0)).is_equal_to(32.0);
        assert_that!(era_dat(2008, 1, 17, 0.0)).is_equal_to(33.0);
        assert_that!(era_dat(2017, 9, 1, 0.0)).is_equal_to(37.0);
    }

    #[test]
    fn cal_to_jd() {
        assert_that!(era_cal_to_jd(2003, 6, 1)).is_equal_to(52791.0);
    }

    #[test]
    fn jd_to_cal() {
        let (y, m, d, f) = era_jd_to_cal(2400000.5, 50123.9999);
        assert_that!(y).is_equal_to(1996);
        assert_that!(m).is_equal_to(2);
        assert_that!(d).is_equal_to(10);
        assert_that!(f)
            .with_abs_tol(1e-7)
            .is_approx_equal_to(0.9999);
    }

    #[test]
    fn sp00() {
        assert_that!(era_sp00(2400000.5, 52541.0).radians().0)
            .with_abs_tol(1e-12)
            .is_approx_equal_to(-0.6216698469981019309e-11);
    }
}
