use crate::math::matrix::Matrix;
use crate::math::matrix_symmetric::MatrixSym;
use crate::math::scalar::Scalar;
use crate::math::vector::Vector;

pub trait Dot<T: Scalar> {
    type Output;

    fn dot(&self, rhs: &Vector<T>) -> Self::Output;
}

pub trait Norm {
    type Output;

    fn norm(&self) -> Self::Output;
}

/// todo: use AddAssign (first implement AddAssign)
impl<T: Scalar> From<MatrixSym<T>> for Matrix<T> {
    fn from(sym: MatrixSym<T>) -> Self {
        let mut mat = Self::zero(sym.rows(), sym.cols());
        for i in 0..sym.rows() {
            for j in 0..=i {
                mat[(i, j)] = sym[(i, j)];
                mat[(j, i)] = sym[(i, j)];
            }
        }
        mat
    }
}

impl<T: Scalar> From<&MatrixSym<T>> for Matrix<T> {
    fn from(sym: &MatrixSym<T>) -> Self {
        let mut mat = Self::zero(sym.rows(), sym.cols());
        for i in 0..sym.rows() {
            for j in 0..=i {
                mat[(i, j)] = sym[(i, j)];
                mat[(j, i)] = sym[(i, j)];
            }
        }
        mat
    }
}
