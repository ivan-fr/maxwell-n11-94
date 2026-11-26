//! Exact outward-rounded arbitrary-precision dyadic interval arithmetic
//! and Pandrosion V31 two-sided order-4 projective radical enclosures.
//!
//! Zero floating-point (`f64`) or `libm` calls are used in any certified bound.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;

/// Exact dyadic rational `m * 2^e`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Dyadic {
    pub m: BigInt,
    pub e: i32,
}

impl Dyadic {
    #[inline]
    pub fn zero() -> Self {
        Self { m: BigInt::zero(), e: 0 }
    }

    #[inline]
    pub fn one() -> Self {
        Self { m: BigInt::one(), e: 0 }
    }

    #[inline]
    pub fn from_i64(v: i64) -> Self {
        if v == 0 { Self::zero() } else { Self { m: BigInt::from(v), e: 0 } }
    }

    #[inline]
    pub fn from_pow2(e: i32) -> Self {
        Self { m: BigInt::one(), e }
    }

    #[inline]
    pub fn is_zero(&self) -> bool {
        self.m.is_zero()
    }

    #[inline]
    pub fn is_positive(&self) -> bool {
        self.m.is_positive()
    }

    #[inline]
    pub fn is_negative(&self) -> bool {
        self.m.is_negative()
    }

    #[inline]
    pub fn neg(&self) -> Self {
        Self { m: -&self.m, e: self.e }
    }

    #[inline]
    pub fn abs(&self) -> Self {
        Self { m: self.m.abs(), e: self.e }
    }

    #[inline]
    pub fn mul_pow2(&self, k: i32) -> Self {
        if self.is_zero() { Self::zero() } else { Self { m: self.m.clone(), e: self.e + k } }
    }

    /// Strip trailing zero bits from `m` (exact normalization).
    pub fn trim(mut self) -> Self {
        if self.m.is_zero() {
            self.e = 0;
            return self;
        }
        if let Some(tz) = self.m.trailing_zeros() {
            if tz > 0 {
                self.m >>= tz as usize;
                self.e += tz as i32;
            }
        }
        self
    }

    /// Round strictly DOWN (`<= exact`) to at most `prec` significant bits.
    pub fn round_down(m: BigInt, e: i32, prec: u32) -> Self {
        if m.is_zero() {
            return Self::zero();
        }
        let bits = m.bits() as u32;
        if bits <= prec {
            return Self { m, e }.trim();
        }
        let sh = (bits - prec) as usize;
        let d = BigInt::one() << sh;
        let q = m.div_floor(&d);
        Self { m: q, e: e + sh as i32 }.trim()
    }

    /// Round strictly UP (`>= exact`) to at most `prec` significant bits.
    pub fn round_up(m: BigInt, e: i32, prec: u32) -> Self {
        if m.is_zero() {
            return Self::zero();
        }
        let bits = m.bits() as u32;
        if bits <= prec {
            return Self { m, e }.trim();
        }
        let sh = (bits - prec) as usize;
        let d = BigInt::one() << sh;
        let q = m.div_ceil(&d);
        Self { m: q, e: e + sh as i32 }.trim()
    }

    /// Exact comparison `self <=> other`.
    pub fn cmp_exact(&self, other: &Self) -> Ordering {
        if self.is_zero() && other.is_zero() {
            return Ordering::Equal;
        }
        if self.is_zero() {
            return if other.is_positive() { Ordering::Less } else { Ordering::Greater };
        }
        if other.is_zero() {
            return if self.is_positive() { Ordering::Greater } else { Ordering::Less };
        }
        match (self.is_positive(), other.is_positive()) {
            (true, false) => return Ordering::Greater,
            (false, true) => return Ordering::Less,
            _ => {}
        }
        // Quick magnitude check via floor(log2(|x|)) = bits - 1 + e
        let mag_a = (self.m.bits() as i64) - 1 + (self.e as i64);
        let mag_b = (other.m.bits() as i64) - 1 + (other.e as i64);
        if mag_a != mag_b {
            return if self.is_positive() { mag_a.cmp(&mag_b) } else { mag_b.cmp(&mag_a) };
        }
        if self.e >= other.e {
            let lhs = &self.m << ((self.e - other.e) as usize);
            lhs.cmp(&other.m)
        } else {
            let rhs = &other.m << ((other.e - self.e) as usize);
            self.m.cmp(&rhs)
        }
    }

    /// Exact addition (no rounding), then round down/up.
    fn add_raw(a: &Self, b: &Self) -> (BigInt, i32) {
        if a.is_zero() { return (b.m.clone(), b.e); }
        if b.is_zero() { return (a.m.clone(), a.e); }
        if a.e >= b.e {
            let sh = (a.e - b.e) as usize;
            ((&a.m << sh) + &b.m, b.e)
        } else {
            let sh = (b.e - a.e) as usize;
            (&a.m + (&b.m << sh), a.e)
        }
    }

    pub fn add_down(a: &Self, b: &Self, prec: u32) -> Self {
        let (m, e) = Self::add_raw(a, b);
        Self::round_down(m, e, prec)
    }

    pub fn add_up(a: &Self, b: &Self, prec: u32) -> Self {
        let (m, e) = Self::add_raw(a, b);
        Self::round_up(m, e, prec)
    }

    pub fn sub_down(a: &Self, b: &Self, prec: u32) -> Self {
        Self::add_down(a, &b.neg(), prec)
    }

    pub fn sub_up(a: &Self, b: &Self, prec: u32) -> Self {
        Self::add_up(a, &b.neg(), prec)
    }

    pub fn mul_exact(a: &Self, b: &Self) -> Self {
        if a.is_zero() || b.is_zero() {
            Self::zero()
        } else {
            Self { m: &a.m * &b.m, e: a.e + b.e }.trim()
        }
    }

    pub fn mul_down(a: &Self, b: &Self, prec: u32) -> Self {
        if a.is_zero() || b.is_zero() { return Self::zero(); }
        Self::round_down(&a.m * &b.m, a.e + b.e, prec)
    }

    pub fn mul_up(a: &Self, b: &Self, prec: u32) -> Self {
        if a.is_zero() || b.is_zero() { return Self::zero(); }
        Self::round_up(&a.m * &b.m, a.e + b.e, prec)
    }

    /// Outward-rounded division of `a / b` (requires `b != 0`).
    fn div_directed(a: &Self, b: &Self, prec: u32, round_up_dir: bool) -> Self {
        assert!(!b.is_zero(), "division by zero in Dyadic");
        if a.is_zero() { return Self::zero(); }
        let (num_m, den_m) = if b.m.is_positive() {
            (a.m.clone(), b.m.clone())
        } else {
            (-&a.m, -&b.m)
        };
        let shift = (prec as i64) + (den_m.bits() as i64) - (num_m.bits() as i64) + 3;
        let (q, exp) = if shift >= 0 {
            let sh = shift as usize;
            let scaled = num_m << sh;
            let q = if round_up_dir { scaled.div_ceil(&den_m) } else { scaled.div_floor(&den_m) };
            (q, a.e - b.e - shift as i32)
        } else {
            let sh = (-shift) as usize;
            let scaled_den = den_m << sh;
            let q = if round_up_dir { num_m.div_ceil(&scaled_den) } else { num_m.div_floor(&scaled_den) };
            (q, a.e - b.e - shift as i32)
        };
        if round_up_dir {
            Self::round_up(q, exp, prec)
        } else {
            Self::round_down(q, exp, prec)
        }
    }

    pub fn div_down(a: &Self, b: &Self, prec: u32) -> Self {
        Self::div_directed(a, b, prec, false)
    }

    pub fn div_up(a: &Self, b: &Self, prec: u32) -> Self {
        Self::div_directed(a, b, prec, true)
    }

    /// Approximate `f64` for display/logging ONLY (never used in proofs!).
    pub fn to_f64_approx(&self) -> f64 {
        if self.is_zero() { return 0.0; }
        let sign = if self.m.is_negative() { -1.0 } else { 1.0 };
        let abs_m = self.m.abs();
        let bits = abs_m.bits() as i32;
        let top_shift = (bits - 53).max(0) as usize;
        let top: u64 = (&abs_m >> top_shift).try_into().unwrap_or(0);
        let eff_e = self.e + top_shift as i32;
        sign * (top as f64) * 2.0f64.powi(eff_e)
    }
}

impl Ord for Dyadic {
    fn cmp(&self, other: &Self) -> Ordering {
        self.cmp_exact(other)
    }
}
impl PartialOrd for Dyadic {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Closed dyadic interval `[lo, hi]` with invariant `lo <= hi`.
#[derive(Clone, Debug)]
pub struct Ival {
    pub lo: Dyadic,
    pub hi: Dyadic,
}

impl Ival {
    pub fn point(x: Dyadic) -> Self {
        Self { lo: x.clone(), hi: x }
    }

    pub fn zero() -> Self {
        Self::point(Dyadic::zero())
    }

    pub fn new(lo: Dyadic, hi: Dyadic) -> Self {
        assert!(lo <= hi, "invalid interval: lo > hi");
        Self { lo, hi }
    }

    /// Exact rational `num / 10^pow10` enclosed in a `prec`-bit dyadic interval.
    pub fn from_decimal_scaled(num: &BigInt, pow10: u32, prec: u32) -> Self {
        let den = BigInt::from(10u32).pow(pow10);
        let a = Dyadic { m: num.clone(), e: 0 };
        let b = Dyadic { m: den, e: 0 };
        Self::new(Dyadic::div_down(&a, &b, prec), Dyadic::div_up(&a, &b, prec))
    }

    /// Parse a decimal string (possibly scientific notation like `-1.23e-12`) to exact dyadic center!
    pub fn parse_decimal_to_dyadic(s: &str, prec: u32) -> Dyadic {
        let s = s.trim();
        let (neg, rest) = if let Some(r) = s.strip_prefix('-') {
            (true, r)
        } else if let Some(r) = s.strip_prefix('+') {
            (false, r)
        } else {
            (false, s)
        };
        let (mant_str, exp10): (&str, i32) = if let Some(idx) = rest.find(|c| c == 'e' || c == 'E') {
            (&rest[..idx], rest[idx + 1..].parse().expect("valid exp"))
        } else {
            (rest, 0)
        };
        let (digits, frac_len) = if let Some(dot) = mant_str.find('.') {
            let int_p = &mant_str[..dot];
            let frac_p = &mant_str[dot + 1..];
            (format!("{}{}", int_p, frac_p), frac_p.len() as i32)
        } else {
            (mant_str.to_string(), 0)
        };
        let mut num = BigInt::parse_bytes(digits.as_bytes(), 10).expect("valid digits");
        if neg { num = -num; }
        let total_exp10 = exp10 - frac_len;
        if total_exp10 >= 0 {
            let p10 = BigInt::from(10u32).pow(total_exp10 as u32);
            Dyadic::round_down(num * p10, 0, prec)
        } else {
            let den = BigInt::from(10u32).pow((-total_exp10) as u32);
            Dyadic::div_down(&Dyadic { m: num, e: 0 }, &Dyadic { m: den, e: 0 }, prec)
        }
    }

    pub fn neg(&self) -> Self {
        Self { lo: self.hi.neg(), hi: self.lo.neg() }
    }

    pub fn add(&self, rhs: &Self, prec: u32) -> Self {
        Self {
            lo: Dyadic::add_down(&self.lo, &rhs.lo, prec),
            hi: Dyadic::add_up(&self.hi, &rhs.hi, prec),
        }
    }

    pub fn sub(&self, rhs: &Self, prec: u32) -> Self {
        Self {
            lo: Dyadic::sub_down(&self.lo, &rhs.hi, prec),
            hi: Dyadic::sub_up(&self.hi, &rhs.lo, prec),
        }
    }

    pub fn mul(&self, rhs: &Self, prec: u32) -> Self {
        let p1_lo = Dyadic::mul_down(&self.lo, &rhs.lo, prec);
        let p2_lo = Dyadic::mul_down(&self.lo, &rhs.hi, prec);
        let p3_lo = Dyadic::mul_down(&self.hi, &rhs.lo, prec);
        let p4_lo = Dyadic::mul_down(&self.hi, &rhs.hi, prec);
        let lo = p1_lo.min(p2_lo).min(p3_lo).min(p4_lo);

        let p1_hi = Dyadic::mul_up(&self.lo, &rhs.lo, prec);
        let p2_hi = Dyadic::mul_up(&self.lo, &rhs.hi, prec);
        let p3_hi = Dyadic::mul_up(&self.hi, &rhs.lo, prec);
        let p4_hi = Dyadic::mul_up(&self.hi, &rhs.hi, prec);
        let hi = p1_hi.max(p2_hi).max(p3_hi).max(p4_hi);
        Self { lo, hi }
    }

    pub fn sqr(&self, prec: u32) -> Self {
        if !self.lo.is_negative() {
            Self {
                lo: Dyadic::mul_down(&self.lo, &self.lo, prec),
                hi: Dyadic::mul_up(&self.hi, &self.hi, prec),
            }
        } else if !self.hi.is_positive() {
            Self {
                lo: Dyadic::mul_down(&self.hi, &self.hi, prec),
                hi: Dyadic::mul_up(&self.lo, &self.lo, prec),
            }
        } else {
            let a = Dyadic::mul_up(&self.lo, &self.lo, prec);
            let b = Dyadic::mul_up(&self.hi, &self.hi, prec);
            Self { lo: Dyadic::zero(), hi: a.max(b) }
        }
    }

    pub fn div(&self, rhs: &Self, prec: u32) -> Self {
        assert!(
            rhs.lo.is_positive() || rhs.hi.is_negative(),
            "interval division straddles zero"
        );
        let p1_lo = Dyadic::div_down(&self.lo, &rhs.lo, prec);
        let p2_lo = Dyadic::div_down(&self.lo, &rhs.hi, prec);
        let p3_lo = Dyadic::div_down(&self.hi, &rhs.lo, prec);
        let p4_lo = Dyadic::div_down(&self.hi, &rhs.hi, prec);
        let lo = p1_lo.min(p2_lo).min(p3_lo).min(p4_lo);

        let p1_hi = Dyadic::div_up(&self.lo, &rhs.lo, prec);
        let p2_hi = Dyadic::div_up(&self.lo, &rhs.hi, prec);
        let p3_hi = Dyadic::div_up(&self.hi, &rhs.lo, prec);
        let p4_hi = Dyadic::div_up(&self.hi, &rhs.hi, prec);
        let hi = p1_hi.max(p2_hi).max(p3_hi).max(p4_hi);
        Self { lo, hi }
    }

    pub fn scale_i64(&self, k: i64, prec: u32) -> Self {
        self.mul(&Self::point(Dyadic::from_i64(k)), prec)
    }

    pub fn abs_max(&self) -> Dyadic {
        self.lo.abs().max(self.hi.abs())
    }

    pub fn width(&self, prec: u32) -> Dyadic {
        Dyadic::sub_up(&self.hi, &self.lo, prec)
    }
}

/// Seed approximation of `sqrt(S)` for a positive dyadic `S = m * 2^e` using pure integer `BigInt::sqrt`,
/// refined by pure-rational Pandrosion V31 projective order-4 steps.
fn pandrosion_v31_seed(s: &Dyadic, prec: u32) -> Dyadic {
    assert!(s.is_positive());
    // Compute a coarse 32-bit seed from leading bits of S (pure integer shift)
    let bits = s.m.bits() as i32;
    let total_log2 = (bits - 1) + s.e;
    let half_e = total_log2.div_euclid(2);
    // Start with a 24-bit integer sqrt seed
    let target_shift = 48 - bits + ((total_log2 & 1) as i32);
    let m_shifted = if target_shift >= 0 {
        &s.m << (target_shift as usize)
    } else {
        &s.m >> ((-target_shift) as usize)
    };
    let mut y = Dyadic {
        m: m_shifted.sqrt(),
        e: half_e - 24,
    };
    if y.is_zero() {
        y = Dyadic::from_pow2(half_e);
    }
    // 4 iterations of Pandrosion V31 projective order-4 map:
    //   Phi_proj(S, y) = (S^2 + 6*y^2*S + y^4) / (4*y*(y^2 + S))
    // Order 4: 24 -> 96 -> 384 -> 1536 bits in 3 steps!
    let six = Dyadic::from_i64(6);
    for _ in 0..4 {
        let y2 = Dyadic::mul_down(&y, &y, prec + 16);
        let y4 = Dyadic::mul_down(&y2, &y2, prec + 16);
        let s2 = Dyadic::mul_down(s, s, prec + 16);
        let mid = Dyadic::mul_down(&six, &Dyadic::mul_down(&y2, s, prec + 16), prec + 16);
        let num = Dyadic::add_down(&Dyadic::add_down(&s2, &mid, prec + 16), &y4, prec + 16);
        let den = Dyadic::mul_down(&y, &Dyadic::add_down(&y2, s, prec + 16), prec + 16).mul_pow2(2);
        y = Dyadic::div_down(&num, &den, prec + 16);
    }
    y
}

/// **Pandrosion V31 Projective Order-4 Two-Sided Radical Enclosure** for `sqrt([s.lo, s.hi])`.
///
/// Uses the exact algebraic identities (valid for all `S > 0, y > 0`):
/// ```text
///   lower(S, y) = 4*y*S*(y^2 + S) / (S^2 + 6*y^2*S + y^4)  <=  sqrt(S)
///   upper(S, y) = (S^2 + 6*y^2*S + y^4) / (4*y*(y^2 + S))  >=  sqrt(S)
/// ```
/// evaluated with strictly outward-rounded dyadic arithmetic, followed by an
/// **independent exact integer squaring assertion** `r_lo^2 <= s.lo` and `r_hi^2 >= s.hi`.
pub fn pandrosion_v31_sqrt(s: &Ival, prec: u32) -> Ival {
    assert!(s.lo.is_positive(), "sqrt of non-positive interval");
    let six = Dyadic::from_i64(6);

    // 1. Lower bound on sqrt(s.lo): round numerator DOWN, denominator UP!
    let yl = pandrosion_v31_seed(&s.lo, prec);
    let yl2_d = Dyadic::mul_down(&yl, &yl, prec);
    let yl2_u = Dyadic::mul_up(&yl, &yl, prec);
    let num_l = Dyadic::mul_down(
        &Dyadic::mul_down(&yl, &s.lo, prec),
        &Dyadic::add_down(&yl2_d, &s.lo, prec),
        prec,
    )
    .mul_pow2(2);
    let den_l = Dyadic::add_up(
        &Dyadic::add_up(
            &Dyadic::mul_up(&s.lo, &s.lo, prec),
            &Dyadic::mul_up(&six, &Dyadic::mul_up(&yl2_u, &s.lo, prec), prec),
            prec,
        ),
        &Dyadic::mul_up(&yl2_u, &yl2_u, prec),
        prec,
    );
    let r_lo = Dyadic::div_down(&num_l, &den_l, prec);

    // 2. Upper bound on sqrt(s.hi): round numerator UP, denominator DOWN!
    let yh = pandrosion_v31_seed(&s.hi, prec);
    let yh2_d = Dyadic::mul_down(&yh, &yh, prec);
    let yh2_u = Dyadic::mul_up(&yh, &yh, prec);
    let num_h = Dyadic::add_up(
        &Dyadic::add_up(
            &Dyadic::mul_up(&s.hi, &s.hi, prec),
            &Dyadic::mul_up(&six, &Dyadic::mul_up(&yh2_u, &s.hi, prec), prec),
            prec,
        ),
        &Dyadic::mul_up(&yh2_u, &yh2_u, prec),
        prec,
    );
    let den_h = Dyadic::mul_down(&yh, &Dyadic::add_down(&yh2_d, &s.hi, prec), prec).mul_pow2(2);
    let r_hi = Dyadic::div_up(&num_h, &den_h, prec);

    // 3. Belt-and-suspenders exact integer squaring certificate!
    assert!(
        r_lo.is_positive() && Dyadic::mul_exact(&r_lo, &r_lo) <= s.lo,
        "Pandrosion V31 lower square-root certificate failed"
    );
    assert!(
        Dyadic::mul_exact(&r_hi, &r_hi) >= s.hi,
        "Pandrosion V31 upper square-root certificate failed"
    );
    Ival::new(r_lo, r_hi)
}
