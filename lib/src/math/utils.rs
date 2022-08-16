use crate::math::matrix::Matrix;
use crate::math::scalar::Scalar;
use crate::math::vector::Vector;

pub fn check_vec_vec<T: Scalar>(s: &str, lhs: &Vector<T>, rhs: &Vector<T>) {
    if lhs.size() != rhs.size() {
        panic!("[{}] Vector dimensions are not the same.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

pub fn check_mat_vec<T: Scalar>(s: &str, lhs: &Matrix<T>, rhs: &Vector<T>) {
    if lhs.cols() != rhs.size() {
        panic!("[{}] Incompatible dimensions.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}
