use crate::tools::math::scalar::Scalar;
use std::iter::FromIterator;
use std::ops::{Index, IndexMut, Mul, MulAssign};

pub type FMatrix = BaseMatrix<f64>;
pub type IMatrix = BaseMatrix<i64>;

#[derive(Clone, Debug, Default, PartialEq)]
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

    fn new_with_value(rows: usize, cols: usize, value: T) -> BaseMatrix<T> {
        if rows <= 0 {
            panic!("Cannot allocate BaseMatrix with rows <= 0");
        } else if cols <= 0 {
            panic!("Cannot allocate BaseMatrix with cols <= 0");
        }

        BaseMatrix::<T> {
            rows,
            cols,
            elem: vec![value; rows * cols],
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new() {
        let mat = FMatrix {
            rows: 5,
            cols: 5,
            elem: vec![0.0; 25],
        };

        assert_eq!(mat, FMatrix::new(5, 5));
    }

    #[test]
    fn new_with_value() {
        let mat = FMatrix {
            rows: 5,
            cols: 5,
            elem: vec![5.0; 25],
        };

        assert_eq!(mat, FMatrix::new_with_value(5, 5, 5.0));
    }

    #[test]
    fn init() {
        let mut mat = FMatrix::new(5, 5);
        mat.init(5.0);

        assert_eq!(mat, FMatrix::new_with_value(5, 5, 5.0));
    }

    #[test]
    fn index() {
        let mat = IMatrix {
            rows: 5,
            cols: 5,
            elem: (0..25).collect(),
        };

        assert_eq!(mat[[3, 0]], 15);
    }

    #[test]
    fn transpose() {
        // square matrix
        let mut mat = IMatrix {
            rows: 3,
            cols: 3,
            elem: (0..9).collect(),
        };
        mat.transpose();

        let mat_ref = IMatrix {
            rows: 3,
            cols: 3,
            elem: vec![0, 3, 6, 1, 4, 7, 2, 5, 8],
        };

        assert_eq!(mat, mat_ref);

        // non-square matrix
        let mut mat = IMatrix {
            rows: 4,
            cols: 3,
            elem: (0..12).collect(),
        };
        mat.transpose();

        let mat_ref = IMatrix {
            rows: 3,
            cols: 4,
            elem: vec![0, 3, 6, 9, 1, 4, 7, 10, 2, 5, 8, 11],
        };

        assert_eq!(mat, mat_ref);
    }

    #[test]
    fn matrix_scalar_mul() {
        let mat = FMatrix::new_with_value(5, 5, 1.0);
        let new_mat = mat * 2.0;

        assert_eq!(new_mat, FMatrix::new_with_value(5, 5, 2.0));
    }

    #[test]
    fn mul_assign() {
        let mut mat = FMatrix::new_with_value(5, 5, 1.0);
        mat *= 2.0;

        assert_eq!(mat, FMatrix::new_with_value(5, 5, 2.0));
    }
}
