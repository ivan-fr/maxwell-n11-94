pub mod coulomb;
pub mod dyadic;
pub mod limiting_poly;

use crate::coulomb::{CertifiedRoot, build_exact_charges, certify_critical_point};
use crate::dyadic::{Dyadic, Ival};
use rayon::prelude::*;
use serde::Deserialize;

#[derive(Deserialize)]
pub struct ChargeJson {
    pub x: String,
    pub y: String,
    pub z: String,
    pub q: String,
}

#[derive(Deserialize)]
pub struct CritJson {
    pub kind: String,
    pub idx: u8,
    pub x: String,
    pub y: String,
    pub z: String,
}

#[derive(Deserialize)]
pub struct WitnessJson {
    pub n_charges: usize,
    pub n_critical_points: usize,
    pub c1: usize,
    pub c2: usize,
    pub decimal_scale_pow10: u32,
    pub charges_num: Vec<ChargeJson>,
    pub critical_points_str: Vec<CritJson>,
}

pub struct FullCertificate {
    pub n_charges: usize,
    pub n_certified: usize,
    pub n_inner: usize,
    pub n_outer: usize,
    pub c1: usize,
    pub c2: usize,
    pub max_contraction_norm: f64,
    pub min_pole_dist: f64,
    pub max_box_radius_inner: f64,
    pub max_box_radius_outer: f64,
}

/// Precision in bits for exact outward-rounded dyadic interval arithmetic (512 bits ~ 154 decimal digits).
pub const CERT_PREC_BITS: u32 = 512;

/// Run the full bit-exact 512-bit Pandrosion V31 + 3D Interval Krawczyk certification
/// of all 94 non-degenerate critical points on the explicit exact-rational 11-charge potential.
pub fn run_full_coulomb_certification(witness_json: &str) -> FullCertificate {
    let w: WitnessJson = serde_json::from_str(witness_json).expect("valid witness JSON");
    assert_eq!(w.n_charges, 11);
    assert_eq!(w.charges_num.len(), 11);
    assert_eq!(w.critical_points_str.len(), 94);

    let raw_charges: Vec<(String, String, String, String)> = w
        .charges_num
        .iter()
        .map(|c| (c.x.clone(), c.y.clone(), c.z.clone(), c.q.clone()))
        .collect();
    let charges = build_exact_charges(&raw_charges, w.decimal_scale_pow10, CERT_PREC_BITS);

    // 1. Certify all 11 positive charges are pairwise distinct!
    for i in 0..charges.len() {
        for j in (i + 1)..charges.len() {
            let distinct = (0..3).any(|k| {
                charges[i].pos[k].hi < charges[j].pos[k].lo
                    || charges[j].pos[k].hi < charges[i].pos[k].lo
            });
            assert!(distinct, "Charges {} and {} overlap!", i, j);
        }
    }

    // 2. Certify every single one of the 94 critical points in parallel!
    let certified: Vec<(String, u8, CertifiedRoot)> = w
        .critical_points_str
        .par_iter()
        .map(|cp| {
            let center: [Dyadic; 3] = [
                Ival::parse_decimal_to_dyadic(&cp.x, CERT_PREC_BITS),
                Ival::parse_decimal_to_dyadic(&cp.y, CERT_PREC_BITS),
                Ival::parse_decimal_to_dyadic(&cp.z, CERT_PREC_BITS),
            ];
            let cr = certify_critical_point(&center, &charges, CERT_PREC_BITS);
            assert_eq!(cr.morse_index, cp.idx, "Morse index mismatch!");
            (cp.kind.clone(), cp.idx, cr)
        })
        .collect();

    // 3. Certify ALL 94 boxes are PAIRWISE DISJOINT!
    for i in 0..certified.len() {
        for j in (i + 1)..certified.len() {
            assert!(
                certified[i].2.box3.disjoint_from(&certified[j].2.box3),
                "Certified boxes {} and {} are NOT disjoint!",
                i,
                j
            );
        }
    }

    // 4. Aggregate counts and verify exact Poincare-Hopf balance c2 - c1 == n - 1 == 10!
    let n_inner = certified.iter().filter(|x| x.0 == "inner").count();
    let n_outer = certified.iter().filter(|x| x.0 == "outer").count();
    let c1 = certified.iter().filter(|x| x.2.morse_index == 1).count();
    let c2 = certified.iter().filter(|x| x.2.morse_index == 2).count();
    assert_eq!(c1 + c2, 94);
    assert_eq!(c2 as isize - c1 as isize, (w.n_charges as isize) - 1);

    let max_contraction_norm = certified
        .iter()
        .map(|x| x.2.contraction_inf_norm)
        .fold(0.0f64, f64::max);
    let min_pole_dist = certified
        .iter()
        .map(|x| x.2.min_pole_dist)
        .fold(f64::INFINITY, f64::min);
    let max_box_radius_inner = certified
        .iter()
        .filter(|x| x.0 == "inner")
        .map(|x| x.2.radius.to_f64_approx())
        .fold(0.0f64, f64::max);
    let max_box_radius_outer = certified
        .iter()
        .filter(|x| x.0 == "outer")
        .map(|x| x.2.radius.to_f64_approx())
        .fold(0.0f64, f64::max);

    FullCertificate {
        n_charges: w.n_charges,
        n_certified: certified.len(),
        n_inner,
        n_outer,
        c1,
        c2,
        max_contraction_norm,
        min_pole_dist,
        max_box_radius_inner,
        max_box_radius_outer,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn certify_core_p27_with_pandsysroots() {
        let cert = limiting_poly::certify_core_p27();
        assert!(cert.max_laplacian_coeff < 1e-13, "P27 must be harmonic");
        assert_eq!(cert.certified_total, 27);
        assert!(cert.complete);
        assert_eq!(cert.real_count, 27);
        assert_eq!(cert.c1, 16);
        assert_eq!(cert.c2, 11);
    }

    #[test]
    fn certify_limiting_p71_with_pandsysroots() {
        let cert = limiting_poly::certify_limiting_p71();
        assert!(cert.max_laplacian_coeff < 1e-13, "P71 must be harmonic");
        assert_eq!(cert.certified_total, 125);
        assert!(cert.complete);
        assert_eq!(cert.real_count, 71);
        assert_eq!(cert.c1, 35);
        assert_eq!(cert.c2, 36);
    }

    #[test]
    fn certify_all_94_coulomb_equilibria_with_pandrosion_v31_krawczyk() {
        let json = include_str!("../data/witness_n11_94.json");
        let res = run_full_coulomb_certification(json);
        assert_eq!(res.n_charges, 11);
        assert_eq!(res.n_certified, 94);
        assert_eq!(res.n_inner, 71);
        assert_eq!(res.n_outer, 23);
        assert_eq!(res.c1, 42);
        assert_eq!(res.c2, 52);
        assert!(res.max_contraction_norm < 1.0);
    }
}
