use crate::linear_algebra::scalar::Scalar;

use std::{
    cmp::{max, min},
    fmt::{Display, Formatter},
    fs::File,
    io::prelude::*,
    ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign},
    ops::{Deref, DerefMut},
};

use serde::{de::DeserializeOwned, Deserialize, Serialize};

pub type FMatrixSym = MatrixSym<f64>;
pub type IMatrixSym = MatrixSym<i64>;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct MatrixSym<T: Scalar> {
    pub n: usize,
    data: Vec<T>,
}

// simple getters
impl<T: Scalar> MatrixSym<T> {
    pub fn size(&self) -> usize {
        self.n * self.n
    }

    pub fn num_elements(&self) -> usize {
        (self.n * (self.n + 1)) / 2
    }

    fn len(&self) -> usize {
        self.data.len()
    }

    fn shape(&self) -> (usize, usize) {
        (self.n, self.n)
    }

    pub fn is_ok(&self) -> bool {
        self.num_elements() == self.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// implement `Deref` so, data can be accessed implicitly
impl<T: Scalar> Deref for MatrixSym<T> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

// implement `DerefMut` so, data can be accessed implicitly
impl<T: Scalar> DerefMut for MatrixSym<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

fn check<T: Scalar>(s: &str, lhs: &MatrixSym<T>, rhs: &MatrixSym<T>) {
    if lhs.shape() != rhs.shape() {
        panic!("[{}] MatrixSym dimensions are not the same.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

impl<T: Scalar> MatrixSym<T> {
    pub fn new(n: usize) -> MatrixSym<T> {
        if n <= 0 {
            panic!("Cannot allocate MatrixSym with n <= 0");
        }

        MatrixSym::<T> {
            n,
            data: Vec::with_capacity(n * (n + 1) / 2),
        }
    }

    pub fn new_with_value(n: usize, value: T) -> MatrixSym<T> {
        if n <= 0 {
            panic!("Cannot allocate MatrixSym with rows <= 0");
        }

        MatrixSym::<T> {
            n,
            data: vec![value; n * (n + 1) / 2],
        }
    }

    pub fn new_from_vec(n: usize, v: &[T]) -> Self {
        if (n * (n + 1)) / 2 != v.len() {
            panic!("[new_from_vec] dimensions are not correct.");
        }

        Self {
            n,
            data: v.to_vec(),
        }
    }

    pub fn zero(n: usize) -> MatrixSym<T> {
        Self::new_with_value(n, T::default())
    }

    pub fn init(&mut self, value: T) {
        self.data = vec![value; self.num_elements()];
    }

    // fn del(&mut self) {
    //     self.n = 0;
    //     self.data.clear();
    // }

    pub fn transpose(&mut self) {}
}

impl<T: Scalar> Display for MatrixSym<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        if self.is_empty() {
            return Ok(());
        }

        for i in 0..self.n {
            for j in 0..self.n {
                write!(f, "{:14.9} ", self[(i, j)])?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------
// Index
// ---------------------------------------------------------------------

impl<T: Scalar> Index<(usize, usize)> for MatrixSym<T> {
    type Output = T;

    /// Calculates index position (=triangle number) and return correct dataent
    /// Note: rather slow to calculate the index every time an dataent needs
    /// to be accessed
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let min = min(index.0, index.1);
        let max = max(index.0, index.1);
        let index = max * (max + 1) / 2 + min;

        &self.data[index]
    }
}

impl<T: Scalar> IndexMut<(usize, usize)> for MatrixSym<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let min = min(index.0, index.1);
        let max = max(index.0, index.1);
        let index = max * (max + 1) / 2 + min;

        &mut self.data[index]
    }
}

impl<T: Scalar + Serialize + DeserializeOwned> MatrixSym<T> {
    pub fn store(&self, name: &str) {
        let mut buffer = File::create(name).expect("Unable to create file");
        write!(
            buffer,
            "{}",
            toml::to_string(self).expect("Unable to serialize MatrixSym")
        )
        .expect("Unable to write to file");
    }

    pub fn retrieve(name: &str) -> Self {
        let mut file = File::open(name).expect("Unable to open file for reading");
        let mut buffer = String::new();
        file.read_to_string(&mut buffer)
            .expect("Unable to read file");
        toml::from_str(&buffer).expect("Unable to deserialize MatrixSym")
    }
}

macro_rules! add (
    ($($add:ident, $add_method:ident,
       $assign:ident, $assign_method:ident,
       $t:ty);* $(;)*) => {$(
    // ---------------------------------------------------------------------
    // Add
    // ---------------------------------------------------------------------

    // MatrixSym = MatrixSym + MatrixSym
    impl $add<MatrixSym<$t>> for MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $add_method(self, rhs: MatrixSym<$t>) -> Self::Output {
            Self::$add_method(self, &rhs)
        }
    }

    // MatrixSym = MatrixSym + &MatrixSym
    impl $add<&MatrixSym<$t>> for MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $add_method(self, rhs: &MatrixSym<$t>) -> Self::Output {
            check("$add", &self, rhs);

            let mut mat = self.clone();
            mat += rhs;

            mat
        }
    }

    // MatrixSym = &MatrixSym + MatrixSym
    impl $add<MatrixSym<$t>> for &MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $add_method(self, rhs: MatrixSym<$t>) -> Self::Output {
            Self::$add_method(self, &rhs)
        }
    }

    // MatrixSym = &MatrixSym + &MatrixSym
    impl $add<&MatrixSym<$t>> for &MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $add_method(self, rhs: &MatrixSym<$t>) -> Self::Output {
            check("$add", self, rhs);

            let mut mat = self.clone();
            mat += rhs;

            mat
        }
    }

    // ---------------------------------------------------------------------
    // AddAssign
    // ---------------------------------------------------------------------
    // allows MatrixSym += MatrixSym
    impl $assign<MatrixSym<$t>> for MatrixSym<$t> {
        fn $assign_method(&mut self, rhs: Self) {
            Self::$assign_method(self, &rhs)
        }
    }

    // allows MatrixSym += &MatrixSym
    impl $assign<&MatrixSym<$t>> for MatrixSym<$t> {
        fn $assign_method(&mut self, rhs: &Self) {
            check("$assign", self, rhs);

            for i in 0..self.n {
                for j in 0..=i {
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

macro_rules! mul_MatrixSym_scalar (
    ($($imp:ident, $method:ident, $t:ty);* $(;)*) => {$(
    // MatrixSym = MatrixSym * scalar
    impl $imp<$t> for MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: $t) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // MatrixSym = MatrixSym * &scalar
    impl $imp<&$t> for MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: &$t) -> Self::Output {
            if self.is_empty() {
                panic!("[Mul] MatrixSym is empty.");
            }

            Self::Output {
                data: self.data.into_iter().map(|v| v * rhs).collect(),
                ..self
            }
        }
    }

    // MatrixSym = &MatrixSym * scalar
    impl $imp<$t> for &MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: $t) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // MatrixSym = &MatrixSym * &scalar
    impl $imp<&$t> for &MatrixSym<$t> {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: &$t) -> Self::Output {
            if self.is_empty() {
                panic!("[Mul] MatrixSym is empty.");
            }

            Self::Output {
                data: self.data.iter().map(|v| v * rhs).collect(),
                ..*self
            }
        }
    }

    // MatrixSym = scalar * MatrixSym
    impl $imp<MatrixSym<$t>> for $t {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: MatrixSym<$t>) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // MatrixSym = scalar * &MatrixSym
    impl $imp<&MatrixSym<$t>> for $t {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: &MatrixSym<$t>) -> Self::Output {
            if rhs.is_empty() {
                panic!("[Mul] MatrixSym is empty.");
            }

            Self::Output {
                data: rhs.data.iter().map(|v| v * self).collect(),
                ..*rhs
            }
        }
    }

    // MatrixSym = &scalar * MatrixSym
    impl $imp<MatrixSym<$t>> for &$t {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: MatrixSym<$t>) -> Self::Output {
            Self::mul(self, &rhs)
        }
    }

    // MatrixSym = &scalar * &MatrixSym
    impl $imp<&MatrixSym<$t>> for &$t {
        type Output = MatrixSym<$t>;

        fn $method(self, rhs: &MatrixSym<$t>) -> Self::Output {
            if rhs.is_empty() {
                panic!("[Mul] MatrixSym is empty.");
            }

            Self::Output {
                data: rhs.data.iter().map(|v| v * self).collect(),
                ..*rhs
            }
        }
    }
    )*};
);

mul_MatrixSym_scalar!(
    Mul, mul, f64;
    Mul, mul, i64;
);

// ---------------------------------------------------------------------
// Sub
// ---------------------------------------------------------------------

// MatrixSym<T> = MatrixSym<T> - MatrixSym<T>
impl<T: Scalar + Sub<Output = T>> Sub<MatrixSym<T>> for MatrixSym<T> {
    type Output = MatrixSym<T>;

    fn sub(self, rhs: MatrixSym<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// MatrixSym<T> = MatrixSym<T> - &MatrixSym<T>
impl<T: Scalar + Sub<Output = T>> Sub<&MatrixSym<T>> for MatrixSym<T> {
    type Output = MatrixSym<T>;

    fn sub(self, rhs: &MatrixSym<T>) -> Self::Output {
        check("Sub", &self, rhs);

        let mut mat = self.clone();
        mat -= rhs;

        mat
    }
}

// MatrixSym<T> = &MatrixSym<T> - MatrixSym<T>
impl<T: Scalar + Sub<Output = T>> Sub<MatrixSym<T>> for &MatrixSym<T> {
    type Output = MatrixSym<T>;

    fn sub(self, rhs: MatrixSym<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// MatrixSym<T> = &MatrixSym<T> - &MatrixSym<T>
impl<T: Scalar + Sub<Output = T>> Sub<&MatrixSym<T>> for &MatrixSym<T> {
    type Output = MatrixSym<T>;

    fn sub(self, rhs: &MatrixSym<T>) -> Self::Output {
        check("Sub", self, rhs);

        let mut mat = self.clone();
        mat -= rhs;

        mat
    }
}

// ---------------------------------------------------------------------
// SubAssign
// ---------------------------------------------------------------------

// allows MatrixSym -= MatrixSym
impl<T: Scalar> SubAssign<MatrixSym<T>> for MatrixSym<T> {
    fn sub_assign(&mut self, rhs: MatrixSym<T>) {
        Self::sub_assign(self, &rhs)
    }
}

// allows MatrixSym -= &MatrixSym
impl<T: Scalar> SubAssign<&MatrixSym<T>> for MatrixSym<T> {
    fn sub_assign(&mut self, rhs: &Self) {
        check("Sub", self, rhs);

        for i in 0..self.n {
            for j in 0..=i {
                self[(i, j)] -= rhs[(i, j)];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const N: usize = 5;
    const N2: usize = (N * (N + 1)) / 2;

    // todo: add cases that panic

    #[test]
    fn new() {
        let mat = FMatrixSym {
            n: N,
            data: Vec::with_capacity(N2),
        };

        assert_eq!(mat, FMatrixSym::new(N));
    }

    #[test]
    fn new_with_value() {
        let mat = FMatrixSym {
            n: N,
            data: vec![5.0; N2],
        };

        assert_eq!(mat, FMatrixSym::new_with_value(N, 5.0));
    }

    #[test]
    fn new_from_vec() {
        let mat = FMatrixSym::new_from_vec(2, &[1.5, 1.5, 1.5]);

        assert_eq!(mat, FMatrixSym::new_with_value(2, 1.5));
    }

    #[test]
    fn zero() {
        let mat = FMatrixSym {
            n: N,
            data: vec![0.0; N2],
        };

        assert_eq!(mat, FMatrixSym::zero(N));
    }

    #[test]
    fn init() {
        let mut mat = FMatrixSym::new(N);
        mat.init(5.0);

        assert_eq!(mat, FMatrixSym::new_with_value(N, 5.0));
    }

    #[test]
    fn index() {
        let mat = IMatrixSym {
            n: N,
            data: (0..N2 as i64).collect(),
        };

        assert_eq!(mat[(3, 0)], 6);
    }

    #[test]
    fn add() {
        let c = FMatrixSym::new_with_value(N, 3.0);

        let a = FMatrixSym::new_with_value(N, 1.0);
        let b = FMatrixSym::new_with_value(N, 2.0);
        assert_eq!(&a + &b, c);
        assert_eq!(a + &b, c);

        let a = FMatrixSym::new_with_value(N, 1.0);
        assert_eq!(&a + b, c);

        let a = FMatrixSym::new_with_value(N, 1.0);
        let b = FMatrixSym::new_with_value(N, 2.0);
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_assign() {
        let c = FMatrixSym::new_with_value(N, 3.0);

        let a = FMatrixSym::new_with_value(N, 1.0);
        let mut b = FMatrixSym::new_with_value(N, 2.0);
        b += &a;
        assert_eq!(b, c);

        let mut b = FMatrixSym::new_with_value(N, 2.0);
        b += a;
        assert_eq!(b, c);
    }

    #[test]
    fn mul_scalar() {
        let mat = FMatrixSym::new_with_value(N, 1.0);
        assert_eq!(&2.0 * &mat, FMatrixSym::new_with_value(N, 2.0));
        assert_eq!(2.0 * &mat, FMatrixSym::new_with_value(N, 2.0));
        assert_eq!(&mat * &2.0, FMatrixSym::new_with_value(N, 2.0));
        assert_eq!(&mat * 2.0, FMatrixSym::new_with_value(N, 2.0));

        let mat = FMatrixSym::new_with_value(N, 1.0);
        assert_eq!(&2.0 * mat, FMatrixSym::new_with_value(N, 2.0));

        let mat = FMatrixSym::new_with_value(N, 1.0);
        assert_eq!(2.0 * mat, FMatrixSym::new_with_value(N, 2.0));

        let mat = FMatrixSym::new_with_value(N, 1.0);
        assert_eq!(mat * &2.0, FMatrixSym::new_with_value(N, 2.0));

        let mat = FMatrixSym::new_with_value(N, 1.0);
        assert_eq!(mat * 2.0, FMatrixSym::new_with_value(N, 2.0));
    }

    #[test]
    fn mul_assign_scalar() {
        let mut mat = FMatrixSym::new_with_value(N, 1.0);

        mat *= 2.0;
        assert_eq!(mat, FMatrixSym::new_with_value(N, 2.0));

        mat *= &2.0;
        assert_eq!(mat, FMatrixSym::new_with_value(N, 4.0));
    }

    #[test]
    fn sub() {
        let c = FMatrixSym::new_with_value(N, 1.0);

        let a = FMatrixSym::new_with_value(N, 3.0);
        let b = FMatrixSym::new_with_value(N, 2.0);
        assert_eq!(&a - &b, c);
        assert_eq!(a - &b, c);

        let a = FMatrixSym::new_with_value(N, 3.0);
        assert_eq!(&a - b, c);

        let a = FMatrixSym::new_with_value(N, 3.0);
        let b = FMatrixSym::new_with_value(N, 2.0);
        assert_eq!(a - b, c);
    }

    #[test]
    fn sub_assign() {
        let c = FMatrixSym::new_with_value(N, 1.0);

        let a = FMatrixSym::new_with_value(N, 2.0);
        let mut b = FMatrixSym::new_with_value(N, 3.0);
        b -= &a;
        assert_eq!(b, c);

        let mut b = FMatrixSym::new_with_value(N, 3.0);
        b -= a;
        assert_eq!(b, c);
    }
}
