use crate::{
    gto_integrals::{e::e, r::r},
    linear_algebra::{
        constants::PI,
        functions::{distance, gaussian_product_center},
    },
};

pub fn eri(
    a: &f64,
    ml_a: &[i16; 3],
    a_origin: &[f64; 3],
    b: &f64,
    ml_b: &[i16; 3],
    b_origin: &[f64; 3],
    c: &f64,
    ml_c: &[i16; 3],
    c_origin: &[f64; 3],
    d: &f64,
    ml_d: &[i16; 3],
    d_origin: &[f64; 3],
) -> f64 {
    let l_a = ml_a[0];
    let m_a = ml_a[1];
    let n_a = ml_a[2];
    let l_b = ml_b[0];
    let m_b = ml_b[1];
    let n_b = ml_b[2];
    let l_c = ml_c[0];
    let m_c = ml_c[1];
    let n_c = ml_c[2];
    let l_d = ml_d[0];
    let m_d = ml_d[1];
    let n_d = ml_d[2];

    let p = a + b;
    let q = c + d;
    let alpha = p * q / (p + q);
    let p_origin = gaussian_product_center(a, a_origin, b, b_origin);
    let q_origin = gaussian_product_center(c, c_origin, d, d_origin);
    let pq_center = [
        p_origin[0] - q_origin[0],
        p_origin[1] - q_origin[1],
        p_origin[2] - q_origin[2],
    ];
    let r_pq = distance(&p_origin, &q_origin);

    let mut val = 0.0;
    for t in 0..=(l_a + l_b) {
        for u in 0..=(m_a + m_b) {
            for v in 0..=(n_a + n_b) {
                for tau in 0..=(l_c + l_d) {
                    for nu in 0..=(m_c + m_d) {
                        for phi in 0..=(n_c + n_d) {
                            val += e(l_a, l_b, t, a_origin[0] - b_origin[0], a, b)
                                * e(m_a, m_b, u, a_origin[1] - b_origin[1], a, b)
                                * e(n_a, n_b, v, a_origin[2] - b_origin[2], a, b)
                                * e(l_c, l_d, tau, c_origin[0] - d_origin[0], c, d)
                                * e(m_c, m_d, nu, c_origin[1] - d_origin[1], c, d)
                                * e(n_c, n_d, phi, c_origin[2] - d_origin[2], c, d)
                                * (-1.0_f64).powf(f64::from(tau + nu + phi))
                                * r(t + tau, u 4+ nu, v + phi, 0, alpha, &pq_center, &r_pq);
                        }
                    }
                }
            }
        }
    }

    val * 2.0 * PI.powf(2.5) / (p * q * (p + q).sqrt())
}
