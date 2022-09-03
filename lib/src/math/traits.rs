use crate::math::matrix::{FMatrix, Matrix};
use crate::math::matrix_symmetric::{FMatrixSym, MatrixSym};
use crate::math::scalar::Scalar;
use crate::math::vector::{FVector, Vector};

use lapack::*;

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

pub trait Diagonalize {
    type Output;

    fn diagonalize(&self) -> Self::Output;
}

// // todo: merge FMatrix FMatrixSym impl (macro?)
// impl Diagonalize for FMatrix {
//     type Output = (FVector, FMatrix);
//
//     fn diagonalize(&self) -> Self::Output {
//         if self.cols() != self.rows() {
//             panic!("[diagonalize] trying to diagonalize a non-square matrix.");
//         }
//
//         let n = self.cols();
//         let mut info = 0;
//         let mut eigenvalues = FVector::zero(n);
//         let mut eigenvectors = self.clone();
//
//         // what is the purpose of this? Is this some kind of buffer for
//         // the diagonalization?
//         let lwork = 3 * n - 1;
//         let mut work = FVector::zero(lwork);
//
//         // todo: careful this assumes a symmetric matrix
//         unsafe {
//             dsyev(
//                 b'V',
//                 b'U',
//                 n as i32,
//                 eigenvectors.as_mut_slice(),
//                 n as i32,
//                 eigenvalues.as_mut_slice(),
//                 work.as_mut_slice(),
//                 lwork as i32,
//                 &mut info,
//             );
//         }
//
//         if info != 0 {
//             panic!("Diagonalization failed with error code: {}", info);
//         }
//
//         (eigenvalues, eigenvectors)
//     }
// }

impl Diagonalize for FMatrixSym {
    type Output = (FVector, FMatrix);

    fn diagonalize(&self) -> Self::Output {
        let n = self.cols();
        let mut info = 0;
        let mut eigenvalues = FVector::zero(n);
        let mut eigenvectors = FMatrix::from(self);

        // what is the purpose of this? Is this some kind of buffer for
        // the diagonalization?
        let lwork = 3 * n - 1;
        let mut work = FVector::zero(lwork);

        unsafe {
            dsyev(
                b'V',
                b'U',
                n as i32,
                eigenvectors.as_mut_slice(),
                n as i32,
                eigenvalues.as_mut_slice(),
                work.as_mut_slice(),
                lwork as i32,
                &mut info,
            );
        }

        if info != 0 {
            panic!("Diagonalization failed with error code: {}", info);
        }

        (eigenvalues, eigenvectors)
    }
}
