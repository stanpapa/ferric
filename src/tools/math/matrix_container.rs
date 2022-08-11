use crate::tools::math::matrix::BaseMatrix;
use crate::tools::math::r#trait::MyTrait;
use std::collections::HashMap;
use std::ops::{Index, IndexMut};

pub type FMatrixContainer = MatrixContainer<f64>;

#[derive(Clone, Default)]
pub struct MatrixContainer<T>
where
    T: MyTrait,
{
    matrices: HashMap<[usize; 2], BaseMatrix<T>>,
}

impl<T: MyTrait> MatrixContainer<T> {
    pub fn new() -> MatrixContainer<T> {
        MatrixContainer::<T>::default()
    }

    pub fn set_matrix(&mut self, index: &[usize; 2], mat: &BaseMatrix<T>) {
        self.matrices.insert(*index, mat.clone());
    }

    pub fn print(&self) {
        for (index, mat) in self.matrices.iter() {
            println!("{:?}", index);
            mat.print()
        }
    }
}

impl<T: MyTrait> Index<[usize; 2]> for MatrixContainer<T> {
    type Output = BaseMatrix<T>;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        match self.matrices.get(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}

impl<T: MyTrait> IndexMut<[usize; 2]> for MatrixContainer<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        match self.matrices.get_mut(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}
