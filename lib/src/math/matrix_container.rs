use crate::math::matrix::Matrix;
use crate::math::scalar::Scalar;
use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::ops::{Index, IndexMut};

pub type FMatrixContainer = MatrixContainer<f64>;

#[derive(Clone, Default)]
pub struct MatrixContainer<T: Scalar> {
    matrices: HashMap<[usize; 2], Matrix<T>>,
}

impl<T: Scalar> MatrixContainer<T> {
    pub fn new() -> MatrixContainer<T> {
        MatrixContainer::<T>::default()
    }

    pub fn insert(&mut self, index: [usize; 2], mat: &Matrix<T>) {
        self.matrices.insert(index, mat.clone());
    }
}

impl<T: Scalar> Display for MatrixContainer<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (index, mat) in &self.matrices {
            writeln!(f, "{:?}\n{}", index, mat)?;
        }

        Ok(())
    }
}

impl<T: Scalar> Index<[usize; 2]> for MatrixContainer<T> {
    type Output = Matrix<T>;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        match self.matrices.get(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}

impl<T: Scalar> IndexMut<[usize; 2]> for MatrixContainer<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        match self.matrices.get_mut(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}
