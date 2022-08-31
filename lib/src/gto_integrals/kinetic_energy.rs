use crate::gto_integrals::overlap::overlap;

pub fn kinetic_energy(
    a: &f64,
    ml_a: &[i16; 3],
    a_origin: &[f64; 3],
    b: &f64,
    ml_b: &[i16; 3],
    b_origin: &[f64; 3],
) -> f64 {
    let l_b = ml_b[0];
    let m_b = ml_b[1];
    let n_b = ml_b[2];

    // double term0 = b*(2*(l2+m2+n2)+3)*overlap(a,lmn1,A,b,lmn2,B);
    let term0 = b
        * (2.0 * f64::from(l_b + m_b + n_b) + 3.0)
        * overlap(a, ml_a, a_origin, b, ml_b, b_origin);

    // double term1 = -2*pow(b,2) *
    // ( overlap( a, lmn1, A, b, [l2+2, m2, n2], B ) +
    // overlap( a, lmn1, A, b, [l2, m2+2, n2], B ) +
    // overlap( a, lmn1, A, b, [l2, m2, n2+2], B ));
    let term1 = -2.0
        * b.powf(2.0)
        * (overlap(a, ml_a, a_origin, b, &[l_b + 2, m_b, n_b], b_origin)
            + overlap(a, ml_a, a_origin, b, &[l_b, m_b + 2, n_b], b_origin)
            + overlap(a, ml_a, a_origin, b, &[l_b, m_b, n_b + 2], b_origin));

    // double term2 = -0.5 * (l2*(l2-1) * overlap( a, lmn1, A, b, [l2-2, m2, n2 ], B ) +
    // m2*(m2-1)*overlap( a, lmn1, A, b, [ l2, m2-2, n2 ], B ) +
    // n2*(n2-1)*overlap( a, lmn1, A, b, [ l2, m2, n2-2 ], B ));
    let term2 = -0.5
        * (f64::from(l_b)
            * f64::from(l_b - 1)
            * overlap(a, ml_a, a_origin, b, &[l_b - 2, m_b, n_b], b_origin)
            + overlap(a, ml_a, a_origin, b, &[l_b, m_b - 2, n_b], b_origin)
            + overlap(a, ml_a, a_origin, b, &[l_b, m_b, n_b - 2], b_origin));

    term0 + term1 + term2
}
