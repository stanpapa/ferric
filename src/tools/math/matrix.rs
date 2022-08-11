use crate::tools::math::scalar::Scalar;
use std::iter::FromIterator;
use std::ops::{Index, IndexMut, Mul, MulAssign};

pub type FMatrix = BaseMatrix<f64>;
pub type IMatrix = BaseMatrix<i64>;

#[derive(Clone, Default)]
pub struct BaseMatrix<T: Scalar> {
    rows: usize,
    cols: usize,
    elem: Vec<T>,
}

// simple getters
impl<T: Scalar> BaseMatrix<T> {
    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    pub fn size(&self) -> usize {
        self.rows * self.cols
    }

    pub fn len(&self) -> usize {
        self.elem.len()
    }

    pub fn shape(&self) -> [usize; 2] {
        [self.rows, self.cols]
    }
}

impl<T: Scalar> BaseMatrix<T> {
    pub fn new(rows: usize, cols: usize) -> BaseMatrix<T> {
        if rows <= 0 {
            panic!("Cannot allocate BaseMatrix with rows <= 0");
        } else if cols <= 0 {
            panic!("Cannot allocate BaseMatrix with cols <= 0");
        }

        BaseMatrix::<T> {
            rows,
            cols,
            elem: vec![T::default(); rows * cols],
        }
    }

    pub fn init(&mut self, value: T) {
        self.elem = vec![value; self.size()];
    }

    fn del(&mut self) {
        self.rows = 0;
        self.cols = 0;
        self.elem.clear();
    }

    pub fn print(&self) {
        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{} ", self[[i, j]]);
            }
            println!()
        }
    }

    pub fn transpose(&mut self) {
        if self.len() == 0 {
            return;
        }

        if self.rows == self.cols {
            // square matrix
            let mut x: T;
            for i in 0..self.rows {
                for j in 0..i {
                    x = self[[i, j]];
                    self[[i, j]] = self[[j, i]];
                    self[[j, i]] = x;
                }
            }
        } else {
            let mat = self.clone();
            self.del();
            *self = BaseMatrix::<T>::new(mat.cols, mat.rows);

            for i in 0..mat.rows {
                for j in 0..mat.cols {
                    self[[j, i]] = mat[[i, j]]
                }
            }
        }
    }
}

impl<T: Scalar> Index<[usize; 2]> for BaseMatrix<T> {
    type Output = T;

    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.elem[index[0] * self.cols + index[1]]
    }
}

impl<T: Scalar> IndexMut<[usize; 2]> for BaseMatrix<T> {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.elem[index[0] * self.cols + index[1]]
    }
}

// allows BaseMatrix * scalar
impl<T: Scalar> Mul<T> for BaseMatrix<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseMatrix<T>;

    fn mul(self, scalar: T) -> Self::Output {
        Self {
            elem: self.elem.into_iter().map(|v| v * scalar).collect(),
            ..self
        }
    }
}

// allows BaseMatrix *= scalar
impl<T: Scalar> MulAssign<T> for BaseMatrix<T> {
    fn mul_assign(&mut self, rhs: T) {
        for x in &mut self.elem {
            *x *= rhs;
        }
    }
}

// todo: add tests
