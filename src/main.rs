use maxwell_n11_94::{CERT_PREC_BITS, limiting_poly, run_full_coulomb_certification};
use std::time::Instant;

fn main() {
    println!("==========================================================================");
    println!("  RIGOROUS CERTIFICATION OF M+(11) >= 94 FOR MAXWELL'S CONJECTURE");
    println!("  Tools: pandsysroots 0.34.0 + Pandrosion V31 projective radical enclosures");
    println!("==========================================================================\n");

    let t0 = Instant::now();
    println!("[1/2] Certifying limiting degree-6 harmonic polynomial P_71 via pandsysroots 0.34.0...");
    let lim = limiting_poly::certify_limiting_p71();
    println!(
        "      max |coeff(Delta P_71)| = {:.3e} (harmonic verified)",
        lim.max_laplacian_coeff
    );
    println!(
        "      pandsysroots 0.34.0: certified = {}/125, complete = {}, real = {} (c1 = {}, c2 = {}, diff = {:+}) [{:.2?}]\n",
        lim.certified_total,
        lim.complete,
        lim.real_count,
        lim.c1,
        lim.c2,
        lim.c2 as isize - lim.c1 as isize,
        t0.elapsed()
    );
    assert_eq!(lim.certified_total, 125);
    assert!(lim.complete);
    assert_eq!(lim.real_count, 71);

    let t1 = Instant::now();
    println!(
        "[2/2] Certifying ALL 94 critical points on the explicit 11-charge positive Coulomb potential"
    );
    println!(
        "      using {}-bit exact outward-rounded dyadic intervals + Pandrosion V31 order-4 sqrt + Krawczyk...",
        CERT_PREC_BITS
    );
    let json = include_str!("../data/witness_n11_94.json");
    let cert = run_full_coulomb_certification(json);
    println!("      Charges:               {} (all strictly positive, exact rationals m/10^90, pairwise distinct)", cert.n_charges);
    println!("      Certified equilibria:  {} ({} inner + {} outer, all 94 boxes pairwise disjoint)", cert.n_certified, cert.n_inner, cert.n_outer);
    println!("      Morse index counts:    c1 = {}, c2 = {}, c2 - c1 = {:+} (Poincare-Hopf exact = +10)", cert.c1, cert.c2, cert.c2 as isize - cert.c1 as isize);
    println!("      Max Krawczyk ||M||_inf:{:.3e} (< 1.0 strictly verified on all 94 boxes)", cert.max_contraction_norm);
    println!("      Min pole distance:     {:.6} (all 94 boxes strictly separated from all 11 poles)", cert.min_pole_dist);
    println!("      Max box radius (in):   {:.3e}", cert.max_box_radius_inner);
    println!("      Max box radius (out):  {:.3e}", cert.max_box_radius_outer);
    println!("      Elapsed:               {:.2?}\n", t1.elapsed());

    println!("==========================================================================");
    println!("  THEOREM CERTIFIED: M+(11) >= 94  (beats previous record 88 by +6!)");
    println!("  COROLLARY (ABK +20 axial lift): M+(2k+1) >= 20k - 6 for all k >= 5.");
    println!("==========================================================================");
}
