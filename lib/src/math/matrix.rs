use crate::math::scalar::Scalar;
use std::fmt::{Display, Formatter};
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign};
use std::ops::{Deref, DerefMut};

pub type FMatrix = Matrix<f64>;
pub type IMatrix = Matrix<i64>;

// // todo: is this OK? or should I use some sort of Matrix trait?
// pub struct MatrixSquare<T: Scalar>(Matrix<T>);
// pub type FMatrixSquare = MatrixSquare<f64>;
// pub type IMatrixSquare = MatrixSquare<i64>;

// pub struct MatrixSym<T: Scalar>(MatrixSquare<T>);
// pub struct MatrixSym<T: Scalar>(Matrix<T>);
// pub type FMatrixSym = MatrixSym<f64>;
// pub type IMatrixSym = MatrixSquare<i64>;

/*
 * note: cannot use macros from forward_ref, because Matrix<T> does
 *       implement Copy...
 */

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Matrix<T: Scalar> {
    rows: usize,
    cols: usize,
    elem: Vec<T>,
}

// simple getters
impl<T: Scalar> Matrix<T> {
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

    fn shape(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn is_ok(&self) -> bool {
        self.size() == self.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// implement `Deref` so, elem can be accessed implicitly
// todo: return slice (deal with lifetimes)
impl<T: Scalar> Deref for Matrix<T> {
    type Target = Vec<T>;

    fn deref(&self) -> &Self::Target {
        &self.elem
    }
}

// WARNING: Do not change len of elem. That will break Vector
// implement `DerefMut` so, elem can be accessed implicitly
impl<T: Scalar> DerefMut for Matrix<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.elem
    }
}

fn check<T: Scalar>(s: &str, lhs: &Matrix<T>, rhs: &Matrix<T>) {
    if lhs.shape() != rhs.shape() {
        panic!("[{}] Matrix dimensions are not the same.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

impl<T: Scalar> Matrix<T> {
    pub fn new(rows: usize, cols: usize) -> Matrix<T> {
        if rows <= 0 {
            panic!("Cannot allocate Matrix with rows <= 0");
        } else if cols <= 0 {
            panic!("Cannot allocate Matrix with cols <= 0");
        }

        Matrix::<T> {
            rows,
            cols,
            elem: Vec::with_capacity(rows * cols),
        }
    }

    pub fn new_with_value(rows: usize, cols: usize, value: T) -> Matrix<T> {
        if rows <= 0 {
            panic!("Cannot allocate Matrix with rows <= 0");
        } else if cols <= 0 {
            panic!("Cannot allocate Matrix with cols <= 0");
        }

        Matrix::<T> {
            rows,
            cols,
            elem: vec![value; rows * cols],
        }
    }

    pub fn new_from_vec(rows: usize, cols: usize, v: &[T]) -> Self {
        if rows * cols != v.len() {
            panic!("[new_from_vec] dimensions are not correct.");
        }

        Self {
            rows,
            cols,
            elem: v.to_vec(),
        }
    }

    pub fn zero(rows: usize, cols: usize) -> Matrix<T> {
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
            *self = Matrix::<T>::new(mat.cols, mat.rows);

            for i in 0..self.rows {
                for j in 0..self.cols {
                    self.elem.push(mat[(j, i)]);
                }
            }
        }
    }

    pub fn transposed(&self) -> Self {
        let mut transposed = self.clone();
        transposed.transpose();
        transposed
    }
}

// ---------------------------------------------------------------------
// Index
// ---------------------------------------------------------------------

impl<T: Scalar> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.elem[index.0 * self.cols + index.1]
    }
}

impl<T: Scalar> IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.elem[index.0 * self.cols + index.1]
    }
}

impl<T: Scalar> Display for Matrix<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() {
            return Ok(());
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                write!(f, "{:14.9} ", self[(i, j)])?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

macro_rules! add (
    ($($add:ident, $add_method:ident,
       $assign:ident, $assign_method:ident,
       $t:ty);* $(;)*) => {$(
    // ---------------------------------------------------------------------
    // Add
    // ---------------------------------------------------------------------

    // Matrix = Matrix + Matrix
    impl $add<Matrix<$t>> for Matrix<$t> {
        type Output = Matrix<$t>;

        fn $add_method(self, rhs: Matrix<$t>) -> Self::Output {
            Self::$add_method(self, &rhs)
        }
    }

    // Matrix = Matrix + &Matrix
    impl $add<&Matrix<$t>> for Matrix<$t> {
        type Output = Matrix<$t>;

        fn $add_method(self, rhs: &Matrix<$t>) -> Self::Output {
            check("$add", &self, rhs);

            let mut mat = self.clone();
            mat += rhs;

            mat
        }
    }

    // Matrix = &Matrix + Matrix
    impl $add<Matrix<$t>> for &Matrix<$t> {
        type Output = Matrix<$t>;

        fn $add_method(self, rhs: Matrix<$t>) -> Self::Output {
            Self::$add_method(self, &rhs)
        }
    }

    // Matrix = &Matrix + &Matrix
    impl $add<&Matrix<$t>> for &Matrix<$t> {
        type Output = Matrix<$t>;

        fn $add_method(self, rhs: &Matrix<$t>) -> Self::Output {
            check("$add", self, rhs);

            let mut mat = self.clone();
            mat += rhs;

            mat
        }
    }

    // ---------------------------------------------------------------------
    // AddAssign
    // ---------------------------------------------------------------------
    // allows Matrix += Matrix
    impl $assign<Matrix<$t>> for Matrix<$t> {
        fn $assign_method(&mut self, rhs: Self) {
            Self::$assign_method(self, &rhs)
        }
    }

    // allows Matrix += &Matrix
    impl $assign<&Matrix<$t>> for Matrix<$t> {
        fn $assign_method(&mut self, rhs: &Self) {
            check("$assign", self, rhs);

            for i in 0..self.rows {
                for j in 0..self.cols {
                    self[(i, j)] += rhs[(i, j)];
                }
            }
        }
    }
    )*};
);

add!(
    Add, add, AddAssign, add_assign, f64;
    Add, add, AddAssign, add_assign, i64;
);

// ---------------------------------------------------------------------
// Mul
// ---------------------------------------------------------------------

macro_rules! mul_matrix_scalar (
    ($($imp:ident, $method:ident, $t:ty);* $(;)*) => {$(
    // Matrix = Matrix * scalar
    impl $imp<$t> for Matrix<$t> {
        type Output = Matrix<$t>;

        fn $method(self, rhs: $t) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // Matrix = Matrix * &scalar
    impl $imp<&$t> for Matrix<$t> {
        type Output = Matrix<$t>;

        fn $method(self, rhs: &$t) -> Self::Output {
            if self.is_empty() {
                panic!("[Mul] Matrix is empty.");
            }

            Self::Output {
                elem: self.elem.into_iter().map(|v| v * rhs).collect(),
                ..self
            }
        }
    }

    // Matrix = &Matrix * scalar
    impl $imp<$t> for &Matrix<$t> {
        type Output = Matrix<$t>;

        fn $method(self, rhs: $t) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // Matrix = &Matrix * &scalar
    impl $imp<&$t> for &Matrix<$t> {
        type Output = Matrix<$t>;

        fn $method(self, rhs: &$t) -> Self::Output {
            if self.is_empty() {
                panic!("[Mul] Matrix is empty.");
            }

            Self::Output {
                elem: self.elem.iter().map(|v| v * rhs).collect(),
                ..*self
            }
        }
    }

    // Matrix = scalar * Matrix
    impl $imp<Matrix<$t>> for $t {
        type Output = Matrix<$t>;

        fn $method(self, rhs: Matrix<$t>) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // Matrix = scalar * &Matrix
    impl $imp<&Matrix<$t>> for $t {
        type Output = Matrix<$t>;

        fn $method(self, rhs: &Matrix<$t>) -> Self::Output {
            if rhs.is_empty() {
                panic!("[Mul] Matrix is empty.");
            }

            Self::Output {
                elem: rhs.elem.iter().map(|v| v * self).collect(),
                ..*rhs
            }
        }
    }

    // Matrix = &scalar * Matrix
    impl $imp<Matrix<$t>> for &$t {
        type Output = Matrix<$t>;

        fn $method(self, rhs: Matrix<$t>) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // Matrix = &scalar * &Matrix
    impl $imp<&Matrix<$t>> for &$t {
        type Output = Matrix<$t>;

        fn $method(self, rhs: &Matrix<$t>) -> Self::Output {
            if rhs.is_empty() {
                panic!("[Mul] Matrix is empty.");
            }

            Self::Output {
                elem: rhs.elem.iter().map(|v| v * self).collect(),
                ..*rhs
            }
        }
    }
    )*};
);

mul_matrix_scalar!(
    Mul, mul, f64;
    Mul, mul, i64;
);

// ---------------------------------------------------------------------
// MulAssign
// ---------------------------------------------------------------------

// allows Matrix *= scalar
// impl<T: Scalar> MulAssign<T> for Matrix<T> {
//     fn mul_assign(&mut self, rhs: T) {
//         Self::mul_assign(self, &rhs)
//     }
// }
//
// // allows Matrix *= &scalar
// impl<T: Scalar> MulAssign<&T> for Matrix<T> {
//     fn mul_assign(&mut self, rhs: &T) {
//         if self.is_empty() {
//             panic!("[Mul] Matrix is empty.");
//         }
//
//         for x in &mut self.elem {
//             *x *= *rhs;
//         }
//     }
// }

// ---------------------------------------------------------------------
// Sub
// ---------------------------------------------------------------------

// Matrix<T> = Matrix<T> - Matrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<Matrix<T>> for Matrix<T> {
    type Output = Matrix<T>;

    fn sub(self, rhs: Matrix<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// Matrix<T> = Matrix<T> - &Matrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<&Matrix<T>> for Matrix<T> {
    type Output = Matrix<T>;

    fn sub(self, rhs: &Matrix<T>) -> Self::Output {
        check("Sub", &self, rhs);

        let mut mat = self.clone();
        mat -= rhs;

        mat
    }
}

// Matrix<T> = &Matrix<T> - Matrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<Matrix<T>> for &Matrix<T> {
    type Output = Matrix<T>;

    fn sub(self, rhs: Matrix<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// Matrix<T> = &Matrix<T> - &Matrix<T>
impl<T: Scalar + Sub<Output = T>> Sub<&Matrix<T>> for &Matrix<T> {
    type Output = Matrix<T>;

    fn sub(self, rhs: &Matrix<T>) -> Self::Output {
        check("Sub", self, rhs);

        let mut mat = self.clone();
        mat -= rhs;

        mat
    }
}

// ---------------------------------------------------------------------
// SubAssign
// ---------------------------------------------------------------------

// allows Matrix -= Matrix
impl<T: Scalar> SubAssign<Matrix<T>> for Matrix<T> {
    fn sub_assign(&mut self, rhs: Matrix<T>) {
        Self::sub_assign(self, &rhs)
    }
}

// allows Matrix -= &Matrix
impl<T: Scalar> SubAssign<&Matrix<T>> for Matrix<T> {
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
    fn new_from_vec() {
        let mat = FMatrix::new_from_vec(2, 2, &[1.5, 1.5, 1.5, 1.5]);

        assert_eq!(mat, FMatrix::new_with_value(2, 2, 1.5));
    }

    #[test]
    fn zero() {
        let mat = FMatrix {
            rows: N,
            cols: N,
            elem: vec![0.0; N2],
        };

        assert_eq!(mat, FMatrix::zero(N, N));
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
        assert_eq!(&2.0 * &mat, FMatrix::new_with_value(N, N, 2.0));
        assert_eq!(2.0 * &mat, FMatrix::new_with_value(N, N, 2.0));
        assert_eq!(&mat * &2.0, FMatrix::new_with_value(N, N, 2.0));
        assert_eq!(&mat * 2.0, FMatrix::new_with_value(N, N, 2.0));

        let mat = FMatrix::new_with_value(N, N, 1.0);
        assert_eq!(&2.0 * mat, FMatrix::new_with_value(N, N, 2.0));

        let mat = FMatrix::new_with_value(N, N, 1.0);
        assert_eq!(2.0 * mat, FMatrix::new_with_value(N, N, 2.0));

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
