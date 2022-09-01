use crate::constants::PI;
use crate::geometry::{distance, Atom};
use crate::gto_integrals::e::e;
use crate::gto_integrals::r::r;
use crate::math::functions::gaussian_product_center;

pub fn nuclear_electron_attraction(
    a: &f64,
    ml_a: &[i16; 3],
    a_origin: &[f64; 3],
    b: &f64,
    ml_b: &[i16; 3],
    b_origin: &[f64; 3],
    atoms: &[Atom],
) -> f64 {
    let l_a = ml_a[0];
    let m_a = ml_a[1];
    let n_a = ml_a[2];

    let l_b = ml_b[0];
    let m_b = ml_b[1];
    let n_b = ml_b[2];

    let p = a + b;
    let p_origin = gaussian_product_center(a, a_origin, b, b_origin);

    let mut val = 0.0;

    for t in 0..=(l_a + l_b) {
        for u in 0..=(m_a + m_b) {
            for v in 0..=(n_a + n_b) {
                // eval = E(l1,l2,t,A[0]-B[0],a,b) * E(m1,m2,u,A[1]-B[1],a,b) * E(n1,n2,v,A[2]-B[2],a,b);
                let eval = e(l_a, l_b, t, a_origin[0] - b_origin[0], a, b)
                    * e(m_a, m_b, u, a_origin[1] - b_origin[1], a, b)
                    * e(n_a, n_b, v, a_origin[2] - b_origin[2], a, b);

                let mut rval = 0.0;
                for atom in atoms {
                    let distance = distance(&p_origin, atom.origin());
                    rval += f64::from(atom.z()) * r(t, u, v, 0, p, &p_origin, &distance);
                }

                val += eval * rval;
            }
        }
    }

    -2.0 * PI / p * val
}
