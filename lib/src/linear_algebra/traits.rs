use crate::linear_algebra::{matrix::Matrix, scalar::Scalar, vector::Vector};

pub trait Dot<T> {
    type Output;

    fn dot(&self, rhs: &T) -> Self::Output;
}

pub trait Norm {
    type Output;

    fn norm(&self) -> Self::Output;
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
