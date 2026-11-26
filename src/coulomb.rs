//! Rigorous 3D Coulomb potential gradient, Hessian, and Krawczyk certification
//! over 512-bit outward-rounded dyadic boxes using Pandrosion V31 radical enclosures.

use crate::dyadic::{Dyadic, Ival, pandrosion_v31_sqrt};
use num_bigint::BigInt;
use num_traits::Signed;

#[derive(Clone, Debug)]
pub struct ChargeBox {
    pub pos: [Ival; 3],
    pub q: Ival,
}

#[derive(Clone, Debug)]
pub struct Box3 {
    pub c: [Ival; 3],
}

impl Box3 {
    pub fn around_center(center: &[Dyadic; 3], radius: &Dyadic, prec: u32) -> Self {
        let mut c = [Ival::zero(), Ival::zero(), Ival::zero()];
        for i in 0..3 {
            c[i] = Ival::new(
                Dyadic::sub_down(&center[i], radius, prec),
                Dyadic::add_up(&center[i], radius, prec),
            );
        }
        Self { c }
    }

    pub fn point(center: &[Dyadic; 3]) -> Self {
        Self {
            c: [
                Ival::point(center[0].clone()),
                Ival::point(center[1].clone()),
                Ival::point(center[2].clone()),
            ],
        }
    }

    pub fn strictly_inside(&self, outer: &Box3) -> bool {
        (0..3).all(|i| self.c[i].lo > outer.c[i].lo && self.c[i].hi < outer.c[i].hi)
    }

    pub fn disjoint_from(&self, other: &Box3) -> bool {
        (0..3).any(|i| self.c[i].hi < other.c[i].lo || other.c[i].hi < self.c[i].lo)
    }
}

/// Evaluate exact outward-rounded enclosures of `grad V(B)` and `Hess V(B)`
/// for `V(x) = sum_j q_j / ||x - a_j||`, using Pandrosion V31 order-4 radical enclosures
/// for every distance `r_j = sqrt(||x - a_j||^2)`.
pub fn eval_grad_hess_interval(
    box3: &Box3,
    charges: &[ChargeBox],
    prec: u32,
) -> ([Ival; 3], [[Ival; 3]; 3], Dyadic) {
    let mut grad = [Ival::zero(), Ival::zero(), Ival::zero()];
    let mut hess = [
        [Ival::zero(), Ival::zero(), Ival::zero()],
        [Ival::zero(), Ival::zero(), Ival::zero()],
        [Ival::zero(), Ival::zero(), Ival::zero()],
    ];
    let mut min_pole_dist_lo: Option<Dyadic> = None;

    for ch in charges {
        let d = [
            box3.c[0].sub(&ch.pos[0], prec),
            box3.c[1].sub(&ch.pos[1], prec),
            box3.c[2].sub(&ch.pos[2], prec),
        ];
        let r2 = d[0].sqr(prec).add(&d[1].sqr(prec), prec).add(&d[2].sqr(prec), prec);
        assert!(r2.lo.is_positive(), "Box touches or encloses a charge pole!");
        // Pandrosion V31 projective order-4 two-sided enclosure of r = sqrt(r2)!
        let r = pandrosion_v31_sqrt(&r2, prec);
        min_pole_dist_lo = Some(match min_pole_dist_lo {
            None => r.lo.clone(),
            Some(prev) => prev.min(r.lo.clone()),
        });
        let r3 = r2.mul(&r, prec);
        let r5 = r3.mul(&r2, prec);
        let q_over_r3 = ch.q.div(&r3, prec);
        let three_q_over_r5 = ch.q.scale_i64(3, prec).div(&r5, prec);

        for a in 0..3 {
            let g_term = q_over_r3.mul(&d[a], prec);
            grad[a] = grad[a].sub(&g_term, prec);
            for b in a..3 {
                let mut h_term = three_q_over_r5.mul(&d[a], prec).mul(&d[b], prec);
                if a == b {
                    h_term = h_term.sub(&q_over_r3, prec);
                }
                hess[a][b] = hess[a][b].add(&h_term, prec);
            }
        }
    }
    hess[1][0] = hess[0][1].clone();
    hess[2][0] = hess[0][2].clone();
    hess[2][1] = hess[1][2].clone();
    (grad, hess, min_pole_dist_lo.unwrap())
}

/// Exact interval 3x3 determinant `det(H)`.
pub fn det3_interval(h: &[[Ival; 3]; 3], prec: u32) -> Ival {
    let c00 = h[1][1].mul(&h[2][2], prec).sub(&h[1][2].mul(&h[2][1], prec), prec);
    let c01 = h[1][0].mul(&h[2][2], prec).sub(&h[1][2].mul(&h[2][0], prec), prec);
    let c02 = h[1][0].mul(&h[2][1], prec).sub(&h[1][1].mul(&h[2][0], prec), prec);
    h[0][0]
        .mul(&c00, prec)
        .sub(&h[0][1].mul(&c01, prec), prec)
        .add(&h[0][2].mul(&c02, prec), prec)
}

/// Approximate 3x3 inverse of the midpoint dyadic matrix `M` (used ONLY as the linear preconditioner `Y`
/// in the Krawczyk operator — any inaccuracy in `Y` can only make the Krawczyk check fail, never falsely pass!).
fn approx_inv3(m: &[[Dyadic; 3]; 3], prec: u32) -> [[Dyadic; 3]; 3] {
    let md = |a: &Dyadic, b: &Dyadic| Dyadic::mul_down(a, b, prec);
    let sb = |a: &Dyadic, b: &Dyadic| Dyadic::sub_down(a, b, prec);
    let ad = |a: &Dyadic, b: &Dyadic| Dyadic::add_down(a, b, prec);

    let c00 = sb(&md(&m[1][1], &m[2][2]), &md(&m[1][2], &m[2][1]));
    let c01 = sb(&md(&m[1][2], &m[2][0]), &md(&m[1][0], &m[2][2]));
    let c02 = sb(&md(&m[1][0], &m[2][1]), &md(&m[1][1], &m[2][0]));
    let c10 = sb(&md(&m[0][2], &m[2][1]), &md(&m[0][1], &m[2][2]));
    let c11 = sb(&md(&m[0][0], &m[2][2]), &md(&m[0][2], &m[2][0]));
    let c12 = sb(&md(&m[0][1], &m[2][0]), &md(&m[0][0], &m[2][1]));
    let c20 = sb(&md(&m[0][1], &m[1][2]), &md(&m[0][2], &m[1][1]));
    let c21 = sb(&md(&m[0][2], &m[1][0]), &md(&m[0][0], &m[1][2]));
    let c22 = sb(&md(&m[0][0], &m[1][1]), &md(&m[0][1], &m[1][0]));

    let det = ad(&ad(&md(&m[0][0], &c00), &md(&m[0][1], &c01)), &md(&m[0][2], &c02));
    assert!(!det.is_zero(), "singular Hessian midpoint");
    // Adjugate is transpose of cofactor matrix
    let dv = |x: &Dyadic| Dyadic::div_down(x, &det, prec);
    [
        [dv(&c00), dv(&c10), dv(&c20)],
        [dv(&c01), dv(&c11), dv(&c21)],
        [dv(&c02), dv(&c12), dv(&c22)],
    ]
}

#[derive(Clone, Debug)]
pub struct CertifiedRoot {
    pub box3: Box3,
    pub radius: Dyadic,
    pub morse_index: u8,
    pub contraction_inf_norm: f64,
    pub min_pole_dist: f64,
}

/// Rigorous 3D Interval Krawczyk Certification for a single critical point around `center`.
///
/// Proves:
/// 1. `K(B) subset int(B)` (strict inclusion in all 3 coordinates),
/// 2. `||I - Y * Hess V(B)||_inf < 1` (contraction => uniqueness and nonsingularity on `B`),
/// 3. `0 not in det(Hess V(B))`, yielding exact Morse index `1` (`det < 0`) or `2` (`det > 0`)
///    since `tr(Hess V) = Delta V = 0` identically outside the charge poles.
pub fn certify_critical_point(
    center: &[Dyadic; 3],
    charges: &[ChargeBox],
    prec: u32,
) -> CertifiedRoot {
    // 1. Evaluate point gradient and Hessian at `center` (as thin intervals)
    let pt_box = Box3::point(center);
    let (g0, h0, _) = eval_grad_hess_interval(&pt_box, charges, prec);

    // Midpoint Hessian for preconditioner Y
    let mut h_mid = [
        [Dyadic::zero(), Dyadic::zero(), Dyadic::zero()],
        [Dyadic::zero(), Dyadic::zero(), Dyadic::zero()],
        [Dyadic::zero(), Dyadic::zero(), Dyadic::zero()],
    ];
    for a in 0..3 {
        for b in 0..3 {
            h_mid[a][b] = Dyadic::add_down(&h0[a][b].lo, &h0[a][b].hi, prec).mul_pow2(-1);
        }
    }
    let y = approx_inv3(&h_mid, prec);
    let y_iv: [[Ival; 3]; 3] = [
        [Ival::point(y[0][0].clone()), Ival::point(y[0][1].clone()), Ival::point(y[0][2].clone())],
        [Ival::point(y[1][0].clone()), Ival::point(y[1][1].clone()), Ival::point(y[1][2].clone())],
        [Ival::point(y[2][0].clone()), Ival::point(y[2][1].clone()), Ival::point(y[2][2].clone())],
    ];

    // Estimate Newton defect `|Y * g0|_inf` to choose an automatic tight box radius!
    let mut yg0 = [Ival::zero(), Ival::zero(), Ival::zero()];
    let mut defect = Dyadic::zero();
    for a in 0..3 {
        let mut s = Ival::zero();
        for b in 0..3 {
            s = s.add(&y_iv[a][b].mul(&g0[b], prec), prec);
        }
        let am = s.abs_max();
        if am > defect {
            defect = am;
        }
        yg0[a] = s;
    }

    // Choose box radius = max(16 * defect, 2^-260) (so Krawczyk contracts by huge margin!)
    let floor_rad = Dyadic::from_pow2(-260);
    let base_rad = defect.mul_pow2(4).max(floor_rad);

    // Try a small ladder of radii (almost always succeeds on the very first radius!)
    for &pow_mul in &[0i32, 4, 8, -4, 12, 16] {
        let rad = base_rad.mul_pow2(pow_mul);
        let box3 = Box3::around_center(center, &rad, prec);
        let (_, hb, min_pole) = eval_grad_hess_interval(&box3, charges, prec);

        // Compute M = I - Y * H(B)
        let mut m = [
            [Ival::zero(), Ival::zero(), Ival::zero()],
            [Ival::zero(), Ival::zero(), Ival::zero()],
            [Ival::zero(), Ival::zero(), Ival::zero()],
        ];
        let mut inf_norm = Dyadic::zero();
        for a in 0..3 {
            let mut row_sum = Dyadic::zero();
            for b in 0..3 {
                let mut prod = Ival::zero();
                for k in 0..3 {
                    prod = prod.add(&y_iv[a][k].mul(&hb[k][b], prec), prec);
                }
                let entry = if a == b {
                    Ival::point(Dyadic::one()).sub(&prod, prec)
                } else {
                    prod.neg()
                };
                row_sum = Dyadic::add_up(&row_sum, &entry.abs_max(), prec);
                m[a][b] = entry;
            }
            if row_sum > inf_norm {
                inf_norm = row_sum;
            }
        }

        // Check contraction ||M||_inf < 1!
        if inf_norm >= Dyadic::one() {
            continue;
        }

        // Compute K(B) = center - Y*g(center) + M * (B - center)
        let diff_iv = Ival::new(rad.neg(), rad.clone());
        let mut k_box_c = [Ival::zero(), Ival::zero(), Ival::zero()];
        for a in 0..3 {
            let mut m_term = Ival::zero();
            for b in 0..3 {
                m_term = m_term.add(&m[a][b].mul(&diff_iv, prec), prec);
            }
            k_box_c[a] = Ival::point(center[a].clone())
                .sub(&yg0[a], prec)
                .add(&m_term, prec);
        }
        let k_box = Box3 { c: k_box_c };
        if !k_box.strictly_inside(&box3) {
            continue;
        }

        // Certify strict sign of det(Hess V(B))!
        let det_iv = det3_interval(&hb, prec);
        let morse_index = if det_iv.hi.is_negative() {
            1u8
        } else if det_iv.lo.is_positive() {
            2u8
        } else {
            continue;
        };

        return CertifiedRoot {
            box3,
            radius: rad,
            morse_index,
            contraction_inf_norm: inf_norm.to_f64_approx(),
            min_pole_dist: min_pole.to_f64_approx(),
        };
    }

    panic!("Krawczyk certification failed for critical point!");
}

pub fn build_exact_charges(raw: &[(String, String, String, String)], pow10: u32, prec: u32) -> Vec<ChargeBox> {
    raw.iter()
        .map(|(xs, ys, zs, qs)| {
            let xn = BigInt::parse_bytes(xs.as_bytes(), 10).expect("valid x int");
            let yn = BigInt::parse_bytes(ys.as_bytes(), 10).expect("valid y int");
            let zn = BigInt::parse_bytes(zs.as_bytes(), 10).expect("valid z int");
            let qn = BigInt::parse_bytes(qs.as_bytes(), 10).expect("valid q int");
            assert!(qn.is_positive(), "Every charge MUST be strictly positive!");
            let q = Ival::from_decimal_scaled(&qn, pow10, prec);
            assert!(q.lo.is_positive(), "Certified lower bound on charge must be > 0!");
            ChargeBox {
                pos: [
                    Ival::from_decimal_scaled(&xn, pow10, prec),
                    Ival::from_decimal_scaled(&yn, pow10, prec),
                    Ival::from_decimal_scaled(&zn, pow10, prec),
                ],
                q,
            }
        })
        .collect()
}
