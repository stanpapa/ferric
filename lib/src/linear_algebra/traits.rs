use crate::linear_algebra::matrix::Matrix;
use crate::linear_algebra::matrix_symmetric::MatrixSym;
use crate::linear_algebra::scalar::Scalar;
use crate::linear_algebra::vector::Vector;

pub trait Dot<T> {
    type Output;

    fn dot(&self, rhs: &T) -> Self::Output;
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
        let mut mat = Self::zero(sym.n, sym.n);
        for i in 0..sym.n {
            for j in 0..=i {
                mat[(i, j)] = sym[(i, j)];
                mat[(j, i)] = sym[(i, j)];
            }
        }
        mat
    }
}

// impl<T: Scalar> From<&mut MatrixSym<T>> for Matrix<T> {
//     fn from(sym: &mut MatrixSym<T>) -> Self {
//         let mut mat = Self::zero(sym.rows(), sym.cols());
//         for i in 0..sym.rows() {
//             for j in 0..=i {
//                 mat[(i, j)] = sym[(i, j)];
//                 mat[(j, i)] = sym[(i, j)];
//             }
//         }
//         mat
//     }
// }

impl<T: Scalar> From<Matrix<T>> for MatrixSym<T> {
    fn from(mat: Matrix<T>) -> Self {
        Self::from(&mat)
    }
}

impl<T: Scalar> From<&Matrix<T>> for MatrixSym<T> {
    fn from(mat: &Matrix<T>) -> Self {
        if mat.rows != mat.cols {
            panic!("Trying to put a non-square matrix into a symmetric matrix.");
        }

        let mut sym = Self::zero(mat.rows);
        for i in 0..mat.rows {
            for j in 0..=i {
                sym[(i, j)] = mat[(i, j)];
            }
        }

        sym
    }
}

pub trait Mat: Clone {
    // getters
    // fn rows(&self) -> usize;
    //
    // fn cols(&self) -> usize;
    //
    // fn size(&self) -> usize;
    //
    // fn len(&self) -> usize;
    //
    // fn shape(&self) -> (usize, usize);
    //
    // // checkers
    // fn is_ok(&self) -> bool;
    // fn is_empty(&self) -> bool {
    //     self.len() == 0
    // }
}

impl<T: Scalar> Mat for Matrix<T> {}
impl<T: Scalar> Mat for MatrixSym<T> {}
