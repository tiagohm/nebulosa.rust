use std::{
    cmp::Ordering,
    f64,
    ops::{Add, Div, Mul, Rem, Sub},
};

/// Adds [a] and [b] exactly, returning the result as two 64-bit floats.
pub fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let x = a + b;
    let mut eb = x - a;
    let mut ea = x - eb;
    eb = b - eb;
    ea = a - ea;
    (x, ea + eb)
}

/// Splits [a] in two aligned parts.
pub fn split(a: f64) -> (f64, f64) {
    let c = 134217729.0 * a;
    let abig = c - a;
    let ah = c - abig;
    (ah, a - ah)
}

/// Multiples [a] and [b] exactly, returning the result as two 64-bit floats.
/// The first is the approximate product (with some floating point error)
/// and the second is the error of the product.
pub fn two_product(a: f64, b: f64) -> (f64, f64) {
    let x = a * b;
    let (ah, al) = split(a);
    let (bh, bl) = split(b);
    let mut y = x - ah * bh;
    y -= al * bh;
    y -= ah * bl;
    (x, al * bl - y)
}

pub trait Zero {
    const ZERO: Self;
}

impl Zero for i8 {
    const ZERO: Self = 0;
}

impl Zero for i16 {
    const ZERO: Self = 0;
}

impl Zero for i32 {
    const ZERO: Self = 0;
}

impl Zero for i64 {
    const ZERO: Self = 0;
}

impl Zero for i128 {
    const ZERO: Self = 0;
}

impl Zero for u8 {
    const ZERO: Self = 0;
}

impl Zero for u16 {
    const ZERO: Self = 0;
}

impl Zero for u32 {
    const ZERO: Self = 0;
}

impl Zero for u64 {
    const ZERO: Self = 0;
}

impl Zero for u128 {
    const ZERO: Self = 0;
}

impl Zero for f32 {
    const ZERO: Self = 0.0;
}

impl Zero for f64 {
    const ZERO: Self = 0.0;
}

/// Computes the modulo where the result is always non-negative.
pub fn pmod<T>(num: T, other: T) -> T
where
    T: Rem<Output = T> + Add<Output = T> + PartialOrd + Copy + Zero,
{
    let rem = num % other;

    if rem < T::ZERO {
        rem + other
    } else {
        rem
    }
}

/// Returns a pair containing the quotient and the remainder when [num] is divided by [other].
pub fn divmod<T>(num: T, other: T) -> (T, T)
where
    T: Rem<Output = T> + Add<Output = T> + Div<Output = T> + PartialOrd + Copy + Zero,
{
    (num / other, pmod(num, other))
}

#[cfg(test)]
mod test {
    use super::pmod;
    use assertor::{assert_that, EqualityAssertion};
    use std::f64::consts::{PI, TAU};

    #[test]
    fn pmod_with_angles() {
        assert_that!(pmod(3.16, PI)).is_equal_to(0.018407346410207026);
        assert_that!(pmod(2.45, PI)).is_equal_to(2.45);
        assert_that!(pmod(-0.018407346410207026, PI)).is_equal_to(3.123185307179586);
        assert_that!(pmod(7.97, TAU)).is_equal_to(1.6868146928204135);
        assert_that!(pmod(5.94, TAU)).is_equal_to(5.94);
        assert_that!(pmod(-1.6868146928204135, TAU)).is_equal_to(4.596370614359173);
    }
}
