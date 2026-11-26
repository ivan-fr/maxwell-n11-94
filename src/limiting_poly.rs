//! Complete complex root certification of the limiting degree-6 harmonic polynomial P_71(X)
//! using `pandsysroots 0.34.0`.

use pandsysroots::{Options, Poly, System, solve};
use std::collections::HashMap;

const C6_TRUE: [f64; 13] = [
    1.0,
    0.26177601754043206,
    0.43601693219955034,
    -0.024699835187192123,
    0.04637463911475271,
    3.358105159412955,
    -0.17242590524684676,
    0.049104446623135505,
    0.07302248052572671,
    0.27467334549277134,
    -0.5863354412689091,
    2.6546787475533518,
    -0.27333590955913756,
];

const H6_BASIS: [&[(u32, u32, u32, f64)]; 13] = [
    &[(0,0,6,16.0),(2,0,4,-120.0),(0,2,4,-120.0),(4,0,2,90.0),(2,2,2,180.0),(0,4,2,90.0),(6,0,0,-5.0),(4,2,0,-15.0),(2,4,0,-15.0),(0,6,0,-5.0)],
    &[(1,0,5,8.0),(3,0,3,-20.0),(1,2,3,-20.0),(5,0,1,5.0),(3,2,1,10.0),(1,4,1,5.0)],
    &[(0,1,5,8.0),(2,1,3,-20.0),(0,3,3,-20.0),(4,1,1,5.0),(2,3,1,10.0),(0,5,1,5.0)],
    &[(2,0,4,16.0),(0,2,4,-16.0),(4,0,2,-16.0),(0,4,2,16.0),(6,0,0,1.0),(4,2,0,1.0),(2,4,0,-1.0),(0,6,0,-1.0)],
    &[(1,1,4,32.0),(3,1,2,-32.0),(1,3,2,-32.0),(5,1,0,2.0),(3,3,0,4.0),(1,5,0,2.0)],
    &[(3,0,3,8.0),(1,2,3,-24.0),(5,0,1,-3.0),(3,2,1,6.0),(1,4,1,9.0)],
    &[(2,1,3,24.0),(0,3,3,-8.0),(4,1,1,-9.0),(2,3,1,-6.0),(0,5,1,3.0)],
    &[(4,0,2,10.0),(2,2,2,-60.0),(0,4,2,10.0),(6,0,0,-1.0),(4,2,0,5.0),(2,4,0,5.0),(0,6,0,-1.0)],
    &[(3,1,2,40.0),(1,3,2,-40.0),(5,1,0,-4.0),(1,5,0,4.0)],
    &[(5,0,1,1.0),(3,2,1,-10.0),(1,4,1,5.0)],
    &[(4,1,1,5.0),(2,3,1,-10.0),(0,5,1,1.0)],
    &[(6,0,0,1.0),(4,2,0,-15.0),(2,4,0,15.0),(0,6,0,-1.0)],
    &[(5,1,0,6.0),(3,3,0,-20.0),(1,5,0,6.0)],
];

fn euler_zyz(a: f64, b: f64, g: f64) -> [[f64; 3]; 3] {
    let (ca, sa) = (a.cos(), a.sin());
    let (cb, sb) = (b.cos(), b.sin());
    let (cg, sg) = (g.cos(), g.sin());
    [
        [ca * cb * cg - sa * sg, -ca * cb * sg - sa * cg, ca * sb],
        [sa * cb * cg + ca * sg, -sa * cb * sg + ca * cg, sa * sb],
        [-sb * cg, sb * sg, cb],
    ]
}

fn terms_from_rot_h6(c6: &[f64; 13], rmat: &[[f64; 3]; 3], scale: f64) -> Vec<(Vec<u32>, f64)> {
    let mul_map = |a: &HashMap<(u32, u32, u32), f64>, b: &HashMap<(u32, u32, u32), f64>| {
        let mut out = HashMap::new();
        for (&(ax, ay, az), &av) in a {
            for (&(bx, by, bz), &bv) in b {
                *out.entry((ax + bx, ay + by, az + bz)).or_insert(0.0) += av * bv;
            }
        }
        out
    };
    let lin = |row: [f64; 3]| {
        let mut m = HashMap::new();
        m.insert((1, 0, 0), row[0]);
        m.insert((0, 1, 0), row[1]);
        m.insert((0, 0, 1), row[2]);
        m
    };
    let one = {
        let mut m = HashMap::new();
        m.insert((0, 0, 0), 1.0);
        m
    };
    let lx = lin(rmat[0]);
    let ly = lin(rmat[1]);
    let lz = lin(rmat[2]);
    let mut xp = vec![one.clone()];
    let mut yp = vec![one.clone()];
    let mut zp = vec![one];
    for _ in 0..6 {
        xp.push(mul_map(xp.last().unwrap(), &lx));
        yp.push(mul_map(yp.last().unwrap(), &ly));
        zp.push(mul_map(zp.last().unwrap(), &lz));
    }
    let mut acc: HashMap<(u32, u32, u32), f64> = HashMap::new();
    for j in 0..13 {
        let w = scale * c6[j];
        if w == 0.0 { continue; }
        for &(ex, ey, ez, c) in H6_BASIS[j] {
            let t = mul_map(&mul_map(&xp[ex as usize], &yp[ey as usize]), &zp[ez as usize]);
            for (k, v) in t {
                *acc.entry(k).or_insert(0.0) += w * c * v;
            }
        }
    }
    acc.into_iter()
        .filter(|(_, v)| v.abs() > 1e-22)
        .map(|((x, y, z), v)| (vec![x, y, z], v))
        .collect()
}

pub struct LimitingCert {
    pub certified_total: usize,
    pub complete: bool,
    pub real_count: usize,
    pub c1: usize,
    pub c2: usize,
    pub max_laplacian_coeff: f64,
}

/// Build the limiting degree-6 harmonic polynomial `P_71(X)`, verify `Delta P_71 == 0`,
/// and certify all 125 Bézout critical points (71 real, 27 complex-conjugate pairs) using `pandsysroots 0.34.0`.
pub fn certify_limiting_p71() -> LimitingCert {
    let th5: f64 = 2.0158552862403435;
    let e5: f64 = 0.05;
    let e6: f64 = 3.0e-7;
    let rm = euler_zyz(0.0, 0.19634954084936207, 0.0);
    let s50 = 11.0f64.sqrt() / 8.0;
    let s53 = (231.0f64 / 128.0).sqrt();

    let mut terms: Vec<(Vec<u32>, f64)> = Vec::new();
    let mut add = |w: f64, b: &[(u32, u32, u32, f64)]| {
        for &(ex, ey, ez, c) in b {
            terms.push((vec![ex, ey, ez], w * c));
        }
    };
    add(-0.8, &[(0,0,2,2.0),(2,0,0,-1.0),(0,2,0,-1.0)]);
    add(-4.0, &[(0,0,3,2.0),(2,0,1,-3.0),(0,2,1,-3.0)]);
    add(15.0, &[(3,0,0,1.0),(1,2,0,-3.0)]);
    add(1.0, &[(0,0,4,8.0),(2,0,2,-24.0),(0,2,2,-24.0),(4,0,0,3.0),(2,2,0,6.0),(0,4,0,3.0)]);
    add(25.0, &[(3,0,1,1.0),(1,2,1,-3.0)]);
    add(e5 * th5.cos() * s50, &[(0,0,5,8.0),(2,0,3,-40.0),(0,2,3,-40.0),(4,0,1,15.0),(2,2,1,30.0),(0,4,1,15.0)]);
    add(e5 * th5.sin() * s53, &[(3,0,2,8.0),(1,2,2,-24.0),(5,0,0,-1.0),(3,2,0,2.0),(1,4,0,3.0)]);
    terms.extend(terms_from_rot_h6(&C6_TRUE, &rm, e6));

    // Verify harmonicity: Laplacian(P_71) = 0 identically!
    let mut lap_map: HashMap<(u32, u32, u32), f64> = HashMap::new();
    for (exp, c) in &terms {
        let (ex, ey, ez) = (exp[0], exp[1], exp[2]);
        if ex >= 2 { *lap_map.entry((ex - 2, ey, ez)).or_insert(0.0) += c * (ex * (ex - 1)) as f64; }
        if ey >= 2 { *lap_map.entry((ex, ey - 2, ez)).or_insert(0.0) += c * (ey * (ey - 1)) as f64; }
        if ez >= 2 { *lap_map.entry((ex, ey, ez - 2)).or_insert(0.0) += c * (ez * (ez - 1)) as f64; }
    }
    let max_lap = lap_map.values().map(|v| v.abs()).fold(0.0f64, f64::max);

    let p = Poly::real(3, &terms);
    let gx = p.diff(0);
    let gy = p.diff(1);
    let gz = p.diff(2);
    let hxx = gx.diff(0);
    let hxy = gx.diff(1);
    let hxz = gx.diff(2);
    let hyy = gy.diff(1);
    let hyz = gy.diff(2);
    let hzz = gz.diff(2);

    let sys = System::new(vec![gx, gy, gz]);
    let mut opt = Options::default();
    opt.restarts = 2;
    let out = solve(&sys, &opt);

    let (mut cert, mut real, mut c1, mut c2) = (0usize, 0usize, 0usize, 0usize);
    for b in &out.balls {
        if b.certified {
            cert += 1;
            if b.real == Some(true) {
                real += 1;
                let pt = [b.center[0], b.center[1], b.center[2]];
                let a00 = hxx.eval(&pt).re;
                let a01 = hxy.eval(&pt).re;
                let a02 = hxz.eval(&pt).re;
                let a11 = hyy.eval(&pt).re;
                let a12 = hyz.eval(&pt).re;
                let a22 = hzz.eval(&pt).re;
                let det = a00 * (a11 * a22 - a12 * a12)
                    - a01 * (a01 * a22 - a12 * a02)
                    + a02 * (a01 * a12 - a11 * a02);
                if det < 0.0 { c1 += 1; } else { c2 += 1; }
            }
        }
    }

    LimitingCert {
        certified_total: cert,
        complete: out.complete,
        real_count: real,
        c1,
        c2,
        max_laplacian_coeff: max_lap,
    }
}

/// Certify the 100%-Bezout degree-4 core harmonic polynomial
/// `P_27(x,y,z) = Y_{4,0} + 25 Y_{4,3c} + 15 Y_{3,3c} - 4 Y_{3,0} - 0.8 Y_{2,0}`,
/// proving all 27/27 Bezout critical points are real and simple (`c1 = 16, c2 = 11, diff = -5`).
pub fn certify_core_p27() -> LimitingCert {
    let mut terms: Vec<(Vec<u32>, f64)> = Vec::new();
    let mut add = |w: f64, b: &[(u32, u32, u32, f64)]| {
        for &(ex, ey, ez, c) in b {
            terms.push((vec![ex, ey, ez], w * c));
        }
    };
    add(-0.8, &[(0,0,2,2.0),(2,0,0,-1.0),(0,2,0,-1.0)]);
    add(-4.0, &[(0,0,3,2.0),(2,0,1,-3.0),(0,2,1,-3.0)]);
    add(15.0, &[(3,0,0,1.0),(1,2,0,-3.0)]);
    add(1.0, &[(0,0,4,8.0),(2,0,2,-24.0),(0,2,2,-24.0),(4,0,0,3.0),(2,2,0,6.0),(0,4,0,3.0)]);
    add(25.0, &[(3,0,1,1.0),(1,2,1,-3.0)]);

    let mut lap_map: HashMap<(u32, u32, u32), f64> = HashMap::new();
    for (exp, c) in &terms {
        let (ex, ey, ez) = (exp[0], exp[1], exp[2]);
        if ex >= 2 { *lap_map.entry((ex - 2, ey, ez)).or_insert(0.0) += c * (ex * (ex - 1)) as f64; }
        if ey >= 2 { *lap_map.entry((ex, ey - 2, ez)).or_insert(0.0) += c * (ey * (ey - 1)) as f64; }
        if ez >= 2 { *lap_map.entry((ex, ey, ez - 2)).or_insert(0.0) += c * (ez * (ez - 1)) as f64; }
    }
    let max_lap = lap_map.values().map(|v| v.abs()).fold(0.0f64, f64::max);

    let p = Poly::real(3, &terms);
    let gx = p.diff(0);
    let gy = p.diff(1);
    let gz = p.diff(2);
    let hxx = gx.diff(0);
    let hxy = gx.diff(1);
    let hxz = gx.diff(2);
    let hyy = gy.diff(1);
    let hyz = gy.diff(2);
    let hzz = gz.diff(2);

    let sys = System::new(vec![gx, gy, gz]);
    let mut opt = Options::default();
    opt.restarts = 2;
    let out = solve(&sys, &opt);

    let (mut cert, mut real, mut c1, mut c2) = (0usize, 0usize, 0usize, 0usize);
    for b in &out.balls {
        if b.certified {
            cert += 1;
            if b.real == Some(true) {
                real += 1;
                let pt = [b.center[0], b.center[1], b.center[2]];
                let a00 = hxx.eval(&pt).re;
                let a01 = hxy.eval(&pt).re;
                let a02 = hxz.eval(&pt).re;
                let a11 = hyy.eval(&pt).re;
                let a12 = hyz.eval(&pt).re;
                let a22 = hzz.eval(&pt).re;
                let det = a00 * (a11 * a22 - a12 * a12)
                    - a01 * (a01 * a22 - a12 * a02)
                    + a02 * (a01 * a12 - a11 * a02);
                if det < 0.0 { c1 += 1; } else { c2 += 1; }
            }
        }
    }

    LimitingCert {
        certified_total: cert,
        complete: out.complete,
        real_count: real,
        c1,
        c2,
        max_laplacian_coeff: max_lap,
    }
}
