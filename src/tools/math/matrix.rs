use crate::tools::math::scalar::Scalar;
use std::iter::FromIterator;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

pub type FMatrix = BaseMatrix<f64>;
pub type IMatrix = BaseMatrix<i64>;

// // todo: is this OK? or should I use some sort of Matrix trait?
// pub struct MatrixSquare<T: Scalar>(BaseMatrix<T>);
// pub type FMatrixSquare = MatrixSquare<f64>;
// pub type IMatrixSquare = MatrixSquare<i64>;
//
// pub struct MatrixSym<T: Scalar>(MatrixSquare<T>);
// pub type FMatrixSym = MatrixSquare<f64>;
// pub type IMatrixSym = MatrixSquare<i64>;

/*
 * note: cannot use macros from forward_ref, because BaseMatrix<T> does
 *       implement Copy...
 */

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

    fn len(&self) -> usize {
        self.elem.len()
    }

    pub fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn is_ok(&self) -> bool {
        self.size() == self.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

fn check<T: Scalar>(s: &str, lhs: &BaseMatrix<T>, rhs: &BaseMatrix<T>) {
    if lhs.shape() != rhs.shape() {
        panic!("[{}] Matrix dimensions are not the same.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
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
            elem: Vec::with_capacity(rows * cols),
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

    pub fn zeros(rows: usize, cols: usize) -> BaseMatrix<T> {
        Self::new_with_value(rows, cols, T::default())
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
        if self.is_empty() {
            println!("Empty matrix. Nothing to print.");
            return;
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{:12.8} ", self[(i, j)]);
            }
            println!()
        }
    }

    pub fn transpose(&mut self) {
        // vector has not been initialized, so we just swap the rows and cols
        if self.is_empty() {
            let rows = self.rows;
            self.rows = self.cols;
            self.cols = rows;

            return;
        }

        if !self.is_ok() {
            panic!("[transpose] Size does not match len!");
        }

        if self.rows == self.cols {
            // square matrix
            let mut x: T;
            for i in 0..self.rows {
                for j in 0..i {
                    x = self[(i, j)];
                    self[(i, j)] = self[(j, i)];
                    self[(j, i)] = x;
                }
            }
        } else {
            let mat = self.clone();
            self.del();
            *self = BaseMatrix::<T>::new(mat.cols, mat.rows);

            for i in 0..self.rows {
                for j in 0..self.cols {
                    self.elem.push(mat[(j, i)]);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------
// Index
// ---------------------------------------------------------------------

impl<T: Scalar> Index<(usize, usize)> for BaseMatrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.elem[index.0 * self.cols + index.1]
    }
}

impl<T: Scalar> IndexMut<(usize, usize)> for BaseMatrix<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.elem[index.0 * self.cols + index.1]
    }
}

// ---------------------------------------------------------------------
// Add
// ---------------------------------------------------------------------

// BaseMatrix<T> = BaseMatrix<T> + BaseMatrix<T>
impl<T: Scalar + Add<Output = T>> Add<BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn add(self, rhs: BaseMatrix<T>) -> Self::Output {
        Self::add(self, &rhs)
    }
}

// BaseMatrix<T> = BaseMatrix<T> + &BaseMatrix<T>
impl<T: Scalar + Add<Output = T>> Add<&BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn add(self, rhs: &BaseMatrix<T>) -> Self::Output {
        check("Add", &self, rhs);

        let mut mat = self.clone();
        mat += rhs;

        mat
    }
}

// BaseMatrix<T> = &BaseMatrix<T> + BaseMatrix<T>
impl<T: Scalar + Add<Output = T>> Add<BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn add(self, rhs: BaseMatrix<T>) -> Self::Output {
        Self::add(self, &rhs)
    }
}

// BaseMatrix<T> = &BaseMatrix<T> + &BaseMatrix<T>
impl<T: Scalar + Add<Output = T>> Add<&BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn add(self, rhs: &BaseMatrix<T>) -> Self::Output {
        check("Add", self, rhs);

        let mut mat = self.clone();
        mat += rhs;

        mat
    }
}

// ---------------------------------------------------------------------
// AddAssign
// ---------------------------------------------------------------------

// allows BaseMatrix += BaseMatrix
impl<T: Scalar> AddAssign<BaseMatrix<T>> for BaseMatrix<T> {
    fn add_assign(&mut self, rhs: Self) {
        Self::add_assign(self, &rhs)
    }
}

// allows BaseMatrix += &BaseMatrix
impl<T: Scalar> AddAssign<&BaseMatrix<T>> for BaseMatrix<T> {
    fn add_assign(&mut self, rhs: &Self) {
        check("AddAssign", self, rhs);

        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] += rhs[(i, j)];
            }
        }
    }
}

// ---------------------------------------------------------------------
// Mul
// ---------------------------------------------------------------------

// BaseMatrix = BaseMatrix * scalar
impl<T: Scalar> Mul<T> for BaseMatrix<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// BaseMatrix = BaseMatrix * &scalar
impl<T: Scalar> Mul<&T> for BaseMatrix<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        if self.is_empty() {
            panic!("[Mul] Matrix is empty.");
        }

        Self::Output {
            elem: self.elem.into_iter().map(|v| v * *rhs).collect(),
            ..self
        }
    }
}

// BaseMatrix = &BaseMatrix * scalar
impl<T: Scalar> Mul<T> for &BaseMatrix<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// BaseMatrix = &BaseMatrix * &scalar
impl<T: Scalar> Mul<&T> for &BaseMatrix<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        if self.is_empty() {
            panic!("[Mul] Matrix is empty.");
        }

        Self::Output {
            elem: self.elem.iter().map(|v| *v * *rhs).collect(),
            ..*self
        }
    }
}

// ---------------------------------------------------------------------
// MulAssign
// ---------------------------------------------------------------------

// allows BaseMatrix *= scalar
impl<T: Scalar> MulAssign<T> for BaseMatrix<T> {
    fn mul_assign(&mut self, rhs: T) {
        Self::mul_assign(self, &rhs)
    }
}

// allows BaseMatrix *= &scalar
impl<T: Scalar> MulAssign<&T> for BaseMatrix<T> {
    fn mul_assign(&mut self, rhs: &T) {
        if self.is_empty() {
            panic!("[Mul] Matrix is empty.");
        }

        for x in &mut self.elem {
            *x *= *rhs;
        }
    }
}

// ---------------------------------------------------------------------
// Sub
// ---------------------------------------------------------------------

// BaseMatrix<T> = BaseMatrix<T> - BaseMatrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn sub(self, rhs: BaseMatrix<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// BaseMatrix<T> = BaseMatrix<T> - &BaseMatrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<&BaseMatrix<T>> for BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn sub(self, rhs: &BaseMatrix<T>) -> Self::Output {
        check("Sub", &self, rhs);

        let mut mat = self.clone();
        mat -= rhs;

        mat
    }
}

// BaseMatrix<T> = &BaseMatrix<T> - BaseMatrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn sub(self, rhs: BaseMatrix<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// BaseMatrix<T> = &BaseMatrix<T> - &BaseMatrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<&BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn sub(self, rhs: &BaseMatrix<T>) -> Self::Output {
        check("Sub", self, rhs);

        let mut mat = self.clone();
        mat -= rhs;

        mat
    }
}

// ---------------------------------------------------------------------
// SubAssign
// ---------------------------------------------------------------------

// allows BaseMatrix -= BaseMatrix
impl<T: Scalar> SubAssign<BaseMatrix<T>> for BaseMatrix<T> {
    fn sub_assign(&mut self, rhs: BaseMatrix<T>) {
        Self::sub_assign(self, &rhs)
    }
}

// allows BaseMatrix -= &BaseMatrix
impl<T: Scalar> SubAssign<&BaseMatrix<T>> for BaseMatrix<T> {
    fn sub_assign(&mut self, rhs: &Self) {
        check("Sub", self, rhs);

        for i in 0..self.rows {
            for j in 0..self.cols {
                self[(i, j)] -= rhs[(i, j)];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const N: usize = 5;
    const N2: usize = N * N;

    // todo: add cases that panic

    #[test]
    fn new() {
        let mat = FMatrix {
            rows: N,
            cols: N,
            elem: Vec::with_capacity(N2),
        };

        assert_eq!(mat, FMatrix::new(N, N));
    }

    #[test]
    fn new_with_value() {
        let mat = FMatrix {
            rows: N,
            cols: N,
            elem: vec![5.0; N2],
        };

        assert_eq!(mat, FMatrix::new_with_value(N, N, 5.0));
    }

    #[test]
    fn zeros() {
        let mat = FMatrix {
            rows: N,
            cols: N,
            elem: vec![0.0; N2],
        };

        assert_eq!(mat, FMatrix::zeros(N, N));
    }

    #[test]
    fn init() {
        let mut mat = FMatrix::new(N, N);
        mat.init(5.0);

        assert_eq!(mat, FMatrix::new_with_value(N, N, 5.0));
    }

    #[test]
    fn index() {
        let mat = IMatrix {
            rows: N,
            cols: N,
            elem: (0..N2 as i64).collect(),
        };

        assert_eq!(mat[(3, 0)], 15);
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
    fn add() {
        let c = FMatrix::new_with_value(N, N, 3.0);

        let a = FMatrix::new_with_value(N, N, 1.0);
        let b = FMatrix::new_with_value(N, N, 2.0);
        assert_eq!(&a + &b, c);
        assert_eq!(a + &b, c);

        let a = FMatrix::new_with_value(N, N, 1.0);
        assert_eq!(&a + b, c);

        let a = FMatrix::new_with_value(N, N, 1.0);
        let b = FMatrix::new_with_value(N, N, 2.0);
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_assign() {
        let c = FMatrix::new_with_value(N, N, 3.0);

        let a = FMatrix::new_with_value(N, N, 1.0);
        let mut b = FMatrix::new_with_value(N, N, 2.0);
        b += &a;
        assert_eq!(b, c);

        let mut b = FMatrix::new_with_value(N, N, 2.0);
        b += a;
        assert_eq!(b, c);
    }

    #[test]
    fn mul_scalar() {
        let mat = FMatrix::new_with_value(N, N, 1.0);
        assert_eq!(&mat * &2.0, FMatrix::new_with_value(N, N, 2.0));

        assert_eq!(&mat * 2.0, FMatrix::new_with_value(N, N, 2.0));

        let mat = FMatrix::new_with_value(N, N, 1.0);
        assert_eq!(mat * &2.0, FMatrix::new_with_value(N, N, 2.0));

        let mat = FMatrix::new_with_value(N, N, 1.0);
        assert_eq!(mat * 2.0, FMatrix::new_with_value(N, N, 2.0));
    }

    #[test]
    fn mul_assign_scalar() {
        let mut mat = FMatrix::new_with_value(N, N, 1.0);

        mat *= 2.0;
        assert_eq!(mat, FMatrix::new_with_value(N, N, 2.0));

        mat *= &2.0;
        assert_eq!(mat, FMatrix::new_with_value(N, N, 4.0));
    }

    #[test]
    fn sub() {
        let c = FMatrix::new_with_value(N, N, 1.0);

        let a = FMatrix::new_with_value(N, N, 3.0);
        let b = FMatrix::new_with_value(N, N, 2.0);
        assert_eq!(&a - &b, c);
        assert_eq!(a - &b, c);

        let a = FMatrix::new_with_value(N, N, 3.0);
        assert_eq!(&a - b, c);

        let a = FMatrix::new_with_value(N, N, 3.0);
        let b = FMatrix::new_with_value(N, N, 2.0);
        assert_eq!(a - b, c);
    }

    #[test]
    fn sub_assign() {
        let c = FMatrix::new_with_value(N, N, 1.0);

        let a = FMatrix::new_with_value(N, N, 2.0);
        let mut b = FMatrix::new_with_value(N, N, 3.0);
        b -= &a;
        assert_eq!(b, c);

        let mut b = FMatrix::new_with_value(N, N, 3.0);
        b -= a;
        assert_eq!(b, c);
    }
}
