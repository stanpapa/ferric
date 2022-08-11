use std::collections::HashMap;
use std::ops::{Index, IndexMut};
use crate::tools::math::matrix::BaseMatrix;

pub struct MatrixContainer<T> {
  matrices: HashMap<[usize; 2], BaseMatrix<T>>,
}

impl<T> Index<[usize; 2]> for MatrixContainer<T> {
    type Output = BaseMatrix<T>;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        match self.matrices.get(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}

impl<T> IndexMut<[usize; 2]> for MatrixContainer<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        match self.matrices.get_mut(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}
