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
        Self::from(&sym)
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

impl<T: Scalar> From<Matrix<T>> for MatrixSym<T> {
    fn from(mat: Matrix<T>) -> Self {
        Self::from(&mat)
    }
}

impl<T: Scalar> From<&Matrix<T>> for MatrixSym<T> {
    fn from(mat: &Matrix<T>) -> Self {
        if mat.rows() != mat.cols() {
            panic!("Trying to put a non-square matrix into a symmetric matrix.");
        }

        let mut sym = Self::new(mat.rows());
        for i in 0..mat.rows() {
            for j in 0..=i {
                sym.push(mat[(i, j)]);
            }
        }

        sym
    }
}

// trait Mat<T: Scalar> {
//     // getters
//     fn rows(&self) -> usize;
//
//     fn cols(&self) -> usize;
//
//     fn size(&self) -> usize;
//
//     fn len(&self) -> usize;
//
//     fn shape(&self) -> (usize, usize);
//
//     // checkers
//     fn is_ok(&self) -> bool;
//     fn is_empty(&self) -> bool {
//         self.len() == 0
//     }
// }
