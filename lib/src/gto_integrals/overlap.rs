use crate::gto_integrals::e::e;
use crate::linear_algebra::constants::PI;

pub fn overlap(
    a: &f64,
    ml_a: &[i16; 3],
    a_origin: &[f64; 3],
    b: &f64,
    ml_b: &[i16; 3],
    b_origin: &[f64; 3],
) -> f64 {
    // double S1 = E(l1,l2,0,A[0]-B[0],a,b);
    // double S2 = E(m1,m2,0,A[1]-B[1],a,b);
    // double S3 = E(n1,n2,0,A[2]-B[2],a,b);
    let s1 = e(ml_a[0], ml_b[0], 0, a_origin[0] - b_origin[0], a, b);
    let s2 = e(ml_a[1], ml_b[1], 0, a_origin[1] - b_origin[1], a, b);
    let s3 = e(ml_a[2], ml_b[2], 0, a_origin[2] - b_origin[2], a, b);

    // return S1*S2*S3*pow(PI/(a+b),1.5);
    s1 * s2 * s3 * (PI / (a + b)).powf(1.5)
}
