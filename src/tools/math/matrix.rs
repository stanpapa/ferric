use std::fmt::Display;
use std::ops::{Index, IndexMut};

pub type FMatrix = BaseMatrix<f64>;
pub type IMatrix = BaseMatrix<i64>;

pub struct BaseMatrix<T> {
    rows: usize,
    cols: usize,
    pub elem: Vec<T>,
}

// simple getters
impl <T> BaseMatrix<T> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn size(&self) -> usize {
        self.rows * self.cols
    }
}

impl<T: Clone + Display> BaseMatrix<T> {
    pub fn new(rows: &usize, cols: &usize) -> BaseMatrix<T> {
        if rows <= &0usize || cols <= &0usize {
            panic!("Cannot allocate BaseMatrix of n <= 0");
        }

        BaseMatrix::<T> {
            rows: *rows,
            cols: *cols,
            elem: Vec::with_capacity((*rows) * (*cols)),
        }
    }

    pub fn init(&mut self, value: &T) {
        self.elem = vec![value.clone(); self.rows * self.cols];
    }

    pub fn print(&self) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{} ", self[[i, j]]);
            }
            println!()
        }
    }
}

impl<T> Index<[usize; 2]> for BaseMatrix<T> {
    type Output = T;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.elem[index[0] * self.rows + index[1]]
    }
}

impl<T> IndexMut<[usize; 2]> for BaseMatrix<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.elem[index[0] * self.rows + index[1]]
    }
}
