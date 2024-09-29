use crate::math::{divmod, two_product, two_sum};
use std::time::{SystemTime, UNIX_EPOCH};

/// Reference epoch (J2000.0), Julian Date.
pub const J2000: f64 = 2451545.0;

/// Reference epoch (B1950.0), Julian Date.
pub const B1950: f64 = 2433282.4235;

/// Seconds per day.
pub const DAYSEC: f64 = 86400.0;

/// Minutes per day.
pub const DAYMIN: f64 = 1440.0;

/// Days per Julian year.
pub const DAYSPERJY: f64 = 365.25;

/// Days per Julian century.
pub const DAYSPERJC: f64 = 36525.0;

/// Days per Julian millennium.
pub const DAYSPERJM: f64 = 365250.0;

/// Julian Date of Modified Julian Date zero.
pub const MJD0: f64 = 2400000.5;

/// Reference epoch (J2000.0), Modified Julian Date.
pub const MJD2000: f64 = 51544.5;

/// 1977 Jan 1.0 as MJD.
pub const MJD1977: f64 = 43144.0;

/// Length of tropical year B1900 (days).
pub const DTY: f64 = 365.242198781;

/// TT minus TAI (s).
pub const TTMINUSTAI: f64 = 32.184;

/// L_G = 1 - d(TT)/d(TCG).
pub const ELG: f64 = 6.969290134E-10;

/// L_B = 1 - d(TDB)/d(TCB).
pub const ELB: f64 = 1.550519768E-8;

/// TDB (s) at TAI 1977/1/1.0.
pub const TDB0: f64 = -6.55E-5;

#[derive(Default, Clone, Copy)]
pub struct Time {
    whole: f64,
    fraction: f64,
}

#[repr(u32)]
pub enum JulianCalendarCutOff {
    None,
    GregorianStart = 2299161,
    GregorianStartEngland = 2361222,
}

trait TimeDelta {
    fn delta(time: &Time) -> f64;
}

impl Time {
    pub fn new(whole: f64, fraction: f64) -> Self {
        Self { whole, fraction }
    }

    pub fn normalized(whole: f64, fraction: f64) -> Self {
        let (whole, fraction) = normalize(whole, fraction, 0.0);
        Self { whole, fraction }
    }

    pub fn from_epoch(epoch: f64, unit: f64, mut whole: f64, mut fraction: f64) -> Self {
        let (a, b) = normalize(epoch, 0.0, unit);
        whole += a;
        fraction += b;

        let extra = f64::round(fraction);
        whole += extra;
        fraction -= extra;

        Self { whole, fraction }
    }

    fn value(&self) -> f64 {
        return self.whole + self.fraction;
    }
}

#[macro_export]
macro_rules! time {
    ($whole: expr, $fraction: expr) => {
        $crate::time::Time::normalized($whole, $fraction)
    };
    ($jd: expr) => {
        $crate::time::Time::normalized($jd, 0.0)
    };
}

pub use time;

/// Unix seconds from 1970-01-01 00:00:00 UTC, ignoring leap seconds.
pub fn time_unix(seconds: f64) -> Time {
    Time::from_epoch(seconds, DAYSEC, 2440588.0, -0.5)
}

pub fn now() -> Time {
    time_unix(
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_secs_f64(),
    )
}

/// Modified Julian Date time format.
pub fn time_mjd(time: f64) -> Time {
    time!(time + MJD0)
}

/// A [Time] represented as [year], [month] and [day].
#[inline(always)]
pub fn time_ymd(year: u32, month: u32, day: u32) -> Time {
    time_ymdhms_cutoff(year, month, day, 0, 0, 0.0, JulianCalendarCutOff::None)
}

/// A [Time] represented as [year], [month] and [day] with [cutoff].
#[inline(always)]
pub fn time_ymd_cutoff(year: u32, month: u32, day: u32, cutoff: JulianCalendarCutOff) -> Time {
    time_ymdhms_cutoff(year, month, day, 0, 0, 0.0, cutoff)
}

/// A [Time] represented as [year], [month], [day], [hour], [minute] and [second].
#[inline(always)]
pub fn time_ymdhms(year: u32, month: u32, day: u32, hour: u32, minute: u32, second: f64) -> Time {
    time_ymdhms_cutoff(
        year,
        month,
        day,
        hour,
        minute,
        second,
        JulianCalendarCutOff::None,
    )
}

/// A [Time] represented as [year], [month], [day], [hour], [minute] and [second].
pub fn time_ymdhms_cutoff(
    year: u32,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: f64,
    cutoff: JulianCalendarCutOff,
) -> Time {
    let (mut y, mut m) = divmod(month - 1, 12);

    y += year;
    m += 1;

    let janfeb = if m <= 2 { 1 } else { 0 };
    let g = y + 4716 - janfeb;
    let f = (m + 9) % 12;
    let e = 1461 * g / 4 + day - 1402;
    let j = e + (153 * f + 2) / 5;

    let fraction = (hour as f64 * 3600.0 + minute as f64 * 60.0 + second) / DAYSEC;

    if j >= cutoff as u32 {
        time!(
            ((j as i32) + (38 - ((g + 184) / 100 * 3 / 4) as i32)) as f64,
            -0.5 + fraction
        )
    } else {
        time!(j as f64, fraction)
    }
}

/// Julian epoch year as floating point value like 2000.0.
pub fn time_julian_epoch(epoch: f64) -> Time {
    time!(2451545.0 + (epoch - 2000.0) * DAYSPERJY)
}

/// Besselian epoch year as floating point value like 1950.0.
pub fn time_besselian_epoch(epoch: f64) -> Time {
    time!(MJD0 + 15019.81352 + (epoch - 1900.0) * DTY)
}

/// Returns the sum of [whole] and [fraction] as two 64-bit floats,
/// with the latter guaranteed to be within -0.5 and 0.5 (inclusive on
/// either side, as the integer is rounded to even).
/// The arithmetic is all done with exact floating point operations so no
/// precision is lost to rounding error. It is assumed the sum is less
/// than about 1E16, otherwise the remainder will be greater than 1.0.
pub fn normalize(whole: f64, fraction: f64, divisor: f64) -> (f64, f64) {
    let (mut sum, mut err) = two_sum(whole, fraction);
    let mut day = f64::round(sum);
    let (extra, mut frac) = two_sum(sum, -day);
    frac += extra + err;

    if (divisor != 0.0 && f64::is_finite(divisor)) {
        let q = sum / divisor;
        let (a, b) = two_product(q, divisor);
        let (c, d) = two_sum(sum, -a);
        (sum, err) = two_sum(q, (c + (d + err - b)) / divisor);
    }

    // Our fraction can now have gotten >0.5 or <-0.5, which means we would
    // loose one bit of precision. So, correct for that.
    day += f64::round(frac);
    let (extra, mut frac) = two_sum(sum, -day);
    frac += extra + err;

    (day, frac)
}

#[cfg(test)]
mod test {
    use super::{
        now, time_besselian_epoch, time_julian_epoch, time_mjd, time_unix, time_ymd, time_ymdhms,
        J2000,
    };
    use assertor::*;

    #[test]
    fn whole_and_fraction() {
        let time = time!(2449353.0, 0.623);
        assert_that!(time.whole).is_equal_to(2449354.0);
        assert_that!(time.fraction).is_equal_to(-0.377);
        assert_that!(time.value()).is_equal_to(2449353.623);
    }

    #[test]
    fn jd() {
        let time = time!(2449353.623);
        assert_that!(time.whole).is_equal_to(2449354.0);
        assert_that!(time.fraction)
            .with_abs_tol(1e-3)
            .is_approx_equal_to(-0.377);
        assert_that!(time.value()).is_equal_to(2449353.623);
    }

    #[test]
    fn unix() {
        assert_that!(time_unix(0.0).value()).is_equal_to(2440587.5);
        assert_that!(time_unix(946684800.0).value()).is_equal_to(J2000 - 0.5);
        assert_that!(now().value()).is_greater_than(2460583.1368056);
    }

    #[test]
    fn mjd() {
        assert_that!(time_mjd(51544.0).value()).is_equal_to(J2000 - 0.5);
    }

    #[test]
    fn ymd() {
        assert_that!(time_ymd(1994, 1, 28).value()).is_equal_to(2449380.5);
        assert_that!(time_ymd(2000, 1, 1).value()).is_equal_to(2451544.5);
    }

    #[test]
    fn ymdhms() {
        assert_that!(time_ymdhms(2003, 4, 25, 14, 7, 35.0).value()).is_equal_to(2452755.088599537);
        assert_that!(time_ymdhms(2000, 1, 1, 12, 0, 0.0).value()).is_equal_to(2451545.0);
    }

    #[test]
    fn julian_epoch() {
        assert_that!(time_julian_epoch(2000.0).value()).is_equal_to(2451545.0);
    }

    #[test]
    fn besselian_epoch() {
        assert_that!(time_besselian_epoch(1950.0).value()).is_equal_to(2433282.42345905);
    }
}
