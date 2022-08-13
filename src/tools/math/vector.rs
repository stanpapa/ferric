use crate::tools::math::scalar::Scalar;
use std::iter::FromIterator;
use std::ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

pub type FVector = BaseVector<f64>;
pub type IVector = BaseVector<i64>;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct BaseVector<T: Scalar> {
    n: usize,
    elem: Vec<T>,
}

// simple getter
impl<T: Scalar> BaseVector<T> {
    pub fn size(&self) -> usize {
        self.n
    }

    fn len(&self) -> usize {
        self.elem.len()
    }

    pub fn is_ok(&self) -> bool {
        self.size() == self.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

fn check<T: Scalar>(s: &str, lhs: &BaseVector<T>, rhs: &BaseVector<T>) {
    if lhs.size() != rhs.size() {
        panic!("[{}] Vector dimensions are not the same.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

impl<T: Scalar> BaseVector<T> {
    pub fn new(n: usize) -> Self {
        if n <= 0 {
            panic!("Cannot allocate BaseVector of n <= 0");
        }

        BaseVector::<T> {
            n,
            elem: Vec::with_capacity(n),
        }
    }

    pub fn new_with_value(n: usize, value: T) -> Self {
        if n <= 0 {
            panic!("Cannot allocate BaseVector of n <= 0");
        }

        BaseVector::<T> {
            n,
            elem: vec![value; n],
        }
    }

    pub fn zeros(n: usize) -> Self {
        Self::new_with_value(n, T::default())
    }

    pub fn init(&mut self, value: T) {
        self.elem = vec![value; self.n];
    }

    fn del(&mut self) {
        self.n = 0;
        self.elem.clear();
    }

    pub fn print(&self) {
        if self.is_empty() {
            println!("Empty vector. Nothing to print.");
            return;
        }

        for i in 0..self.n {
            println!("{} ", self[i]);
        }
    }
}

// ---------------------------------------------------------------------
// Index
// ---------------------------------------------------------------------

impl<T: Scalar> Index<usize> for BaseVector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.elem[index]
    }
}

impl<T: Scalar> IndexMut<usize> for BaseVector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.elem[index]
    }
}

// ---------------------------------------------------------------------
// Add
// ---------------------------------------------------------------------

// BaseVector<T> = BaseVector<T> + BaseVector<T>
impl<T: Scalar + Add<Output = T>> Add<BaseVector<T>> for BaseVector<T> {
    type Output = BaseVector<T>;

    fn add(self, rhs: BaseVector<T>) -> Self::Output {
        Self::add(self, &rhs)
    }
}

// BaseVector<T> = BaseVector<T> + &BaseVector<T>
impl<T: Scalar + Add<Output = T>> Add<&BaseVector<T>> for BaseVector<T> {
    type Output = BaseVector<T>;

    fn add(self, rhs: &BaseVector<T>) -> Self::Output {
        check("Add", &self, rhs);

        let mut vec = self.clone();
        vec += rhs;

        vec
    }
}

// BaseVector<T> = &BaseVector<T> + BaseVector<T>
impl<T: Scalar + Add<Output = T>> Add<BaseVector<T>> for &BaseVector<T> {
    type Output = BaseVector<T>;

    fn add(self, rhs: BaseVector<T>) -> Self::Output {
        Self::add(self, &rhs)
    }
}

// BaseVector<T> = &BaseVector<T> + &BaseVector<T>
impl<T: Scalar + Add<Output = T>> Add<&BaseVector<T>> for &BaseVector<T> {
    type Output = BaseVector<T>;

    fn add(self, rhs: &BaseVector<T>) -> Self::Output {
        check("Add", self, rhs);

        let mut vec = self.clone();
        vec += rhs;

        vec
    }
}

// ---------------------------------------------------------------------
// AddAssign
// ---------------------------------------------------------------------

// allows BaseVector += BaseVector
impl<T: Scalar> AddAssign<BaseVector<T>> for BaseVector<T> {
    fn add_assign(&mut self, rhs: Self) {
        Self::add_assign(self, &rhs)
    }
}

// allows BaseVector += &BaseVector
impl<T: Scalar> AddAssign<&BaseVector<T>> for BaseVector<T> {
    fn add_assign(&mut self, rhs: &Self) {
        check("AddAssign", self, rhs);

        for i in 0..self.n {
            self[i] += rhs[i];
        }
    }
}

// ---------------------------------------------------------------------
// Mul
// ---------------------------------------------------------------------

// BaseVector = BaseVector * scalar
impl<T: Scalar> Mul<T> for BaseVector<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseVector<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// BaseVector = BaseVector * &scalar
impl<T: Scalar> Mul<&T> for BaseVector<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseVector<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        Self::Output {
            elem: self.elem.into_iter().map(|v| v * *rhs).collect(),
            ..self
        }
    }
}

// BaseVector = &BaseVector * scalar
impl<T: Scalar> Mul<T> for &BaseVector<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseVector<T>;

    fn mul(self, rhs: T) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// BaseVector = &BaseVector * &scalar
impl<T: Scalar> Mul<&T> for &BaseVector<T>
where
    Vec<T>: FromIterator<<T as Mul>::Output>,
{
    type Output = BaseVector<T>;

    fn mul(self, rhs: &T) -> Self::Output {
        Self::Output {
            elem: self.elem.iter().map(|v| *v * *rhs).collect(),
            ..*self
        }
    }
}

// todo: make it work with generics
// https://users.rust-lang.org/t/implementing-generic-trait-with-local-struct-on-local-trait/23225
// BaseVector = scalar * BaseVector
// impl<T: Scalar> Mul<BaseVector<T>> for T {
//     type Output = BaseVector<T>;
//
//     fn mul(self, rhs: BaseVector<T>) -> Self::Output {
//         Self::Output {
//             elem: rhs.elem.into_iter().map(|v| v * self).collect(),
//             ..rhs
//         }
//     }
// }

// FVector = scalar * FVector
impl Mul<FVector> for f64 {
    type Output = FVector;

    fn mul(self, rhs: FVector) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// FVector = scalar * &FVector
impl Mul<&FVector> for f64 {
    type Output = FVector;

    fn mul(self, rhs: &FVector) -> Self::Output {
        Self::Output {
            elem: rhs.elem.iter().map(|v| v * self).collect(),
            ..*rhs
        }
    }
}

// FVector = &scalar * FVector
impl Mul<FVector> for &f64 {
    type Output = FVector;

    fn mul(self, rhs: FVector) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// FVector = &scalar * &FVector
impl Mul<&FVector> for &f64 {
    type Output = FVector;

    fn mul(self, rhs: &FVector) -> Self::Output {
        Self::Output {
            elem: rhs.elem.iter().map(|v| v * self).collect(),
            ..*rhs
        }
    }
}

// ---------------------------------------------------------------------
// MulAssign
// ---------------------------------------------------------------------

// BaseVector *= scalar
impl<T: Scalar> MulAssign<T> for BaseVector<T> {
    fn mul_assign(&mut self, rhs: T) {
        Self::mul_assign(self, &rhs)
    }
}

// BaseVector *= &scalar
impl<T: Scalar> MulAssign<&T> for BaseVector<T> {
    fn mul_assign(&mut self, rhs: &T) {
        for x in &mut self.elem {
            *x *= *rhs;
        }
    }
}

// ---------------------------------------------------------------------
// Sub
// ---------------------------------------------------------------------

// BaseVector<T> = BaseVector<T> - BaseVector<T>
impl<T: Scalar + Sub<Output = T>> Sub<BaseVector<T>> for BaseVector<T> {
    type Output = BaseVector<T>;

    fn sub(self, rhs: BaseVector<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// BaseVector<T> = BaseVector<T> - &BaseVector<T>
impl<T: Scalar + Sub<Output = T>> Sub<&BaseVector<T>> for BaseVector<T> {
    type Output = BaseVector<T>;

    fn sub(self, rhs: &BaseVector<T>) -> Self::Output {
        check("Sub", &self, rhs);

        let mut vec = self.clone();
        vec -= rhs;

        vec
    }
}

// BaseVector<T> = &BaseVector<T> - BaseVector<T>
impl<T: Scalar + Sub<Output = T>> Sub<BaseVector<T>> for &BaseVector<T> {
    type Output = BaseVector<T>;

    fn sub(self, rhs: BaseVector<T>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// BaseVector<T> = &BaseVector<T> - &BaseVector<T>
impl<T: Scalar + Sub<Output = T>> Sub<&BaseVector<T>> for &BaseVector<T> {
    type Output = BaseVector<T>;

    fn sub(self, rhs: &BaseVector<T>) -> Self::Output {
        check("Sub", self, rhs);

        let mut vec = self.clone();
        vec -= rhs;

        vec
    }
}

// ---------------------------------------------------------------------
// SubAssign
// ---------------------------------------------------------------------

// allows BaseVector -= BaseVector
impl<T: Scalar> SubAssign<BaseVector<T>> for BaseVector<T> {
    fn sub_assign(&mut self, rhs: BaseVector<T>) {
        Self::sub_assign(self, &rhs);
    }
}

// allows BaseVector -= &BaseVector
impl<T: Scalar> SubAssign<&BaseVector<T>> for BaseVector<T> {
    fn sub_assign(&mut self, rhs: &Self) {
        check("Sub", self, rhs);

        for i in 0..self.n {
            self[i] -= rhs[i];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const N: usize = 5;

    #[test]
    fn new() {
        let vec = FVector {
            n: N,
            elem: Vec::with_capacity(N),
        };

        assert_eq!(vec, FVector::new(N));
    }

    #[test]
    fn new_with_value() {
        let vec = FVector {
            n: N,
            elem: vec![5.0; N],
        };

        assert_eq!(vec, FVector::new_with_value(N, 5.0));
    }

    #[test]
    fn zeros() {
        let vec = FVector {
            n: N,
            elem: vec![0.0; N],
        };

        assert_eq!(vec, FVector::zeros(N));
    }

    #[test]
    fn init() {
        let mut vec = FVector::new(N);
        vec.init(5.0);

        assert_eq!(vec, FVector::new_with_value(N, 5.0));
    }

    #[test]
    fn index() {
        let vec = IVector {
            n: N,
            elem: (0..5).collect(),
        };

        assert_eq!(vec[3], 3);
    }

    #[test]
    fn add() {
        let c = FVector::new_with_value(N, 3.0);

        let a = FVector::new_with_value(N, 1.0);
        let b = FVector::new_with_value(N, 2.0);
        assert_eq!(&a + &b, c);
        assert_eq!(a + &b, c);

        let a = FVector::new_with_value(N, 1.0);
        assert_eq!(&a + b, c);

        let a = FVector::new_with_value(N, 1.0);
        let b = FVector::new_with_value(N, 2.0);
        assert_eq!(a + b, c);
    }

    #[test]
    fn add_assign() {
        let c = FVector::new_with_value(N, 3.0);

        let a = FVector::new_with_value(N, 1.0);
        let mut b = FVector::new_with_value(N, 2.0);
        b += &a;
        assert_eq!(b, c);

        let mut b = FVector::new_with_value(N, 2.0);
        b += a;
        assert_eq!(b, c);
    }

    #[test]
    fn mul_scalar() {
        let vec = FVector::new_with_value(N, 1.0);
        assert_eq!(&2.0 * &vec, FVector::new_with_value(N, 2.0));
        assert_eq!(2.0 * &vec, FVector::new_with_value(N, 2.0));
        assert_eq!(&vec * &2.0, FVector::new_with_value(N, 2.0));
        assert_eq!(&vec * 2.0, FVector::new_with_value(N, 2.0));

        let vec = FVector::new_with_value(N, 1.0);
        assert_eq!(&2.0 * vec, FVector::new_with_value(N, 2.0));

        let vec = FVector::new_with_value(N, 1.0);
        assert_eq!(2.0 * vec, FVector::new_with_value(N, 2.0));

        let vec = FVector::new_with_value(N, 1.0);
        assert_eq!(vec * &2.0, FVector::new_with_value(N, 2.0));

        let vec = FVector::new_with_value(N, 1.0);
        assert_eq!(vec * 2.0, FVector::new_with_value(N, 2.0));
    }

    #[test]
    fn mul_assign_scalar() {
        let mut vec = FVector::new_with_value(N, 1.0);

        vec *= 2.0;
        assert_eq!(vec, FVector::new_with_value(N, 2.0));

        vec *= &2.0;
        assert_eq!(vec, FVector::new_with_value(N, 4.0));
    }

    #[test]
    fn sub() {
        let c = FVector::new_with_value(N, 1.0);

        let a = FVector::new_with_value(N, 3.0);
        let b = FVector::new_with_value(N, 2.0);
        assert_eq!(&a - &b, c);
        assert_eq!(a - &b, c);

        let a = FVector::new_with_value(N, 3.0);
        assert_eq!(&a - b, c);

        let a = FVector::new_with_value(N, 3.0);
        let b = FVector::new_with_value(N, 2.0);
        assert_eq!(a - b, c);
    }

    #[test]
    fn sub_assign() {
        let c = FVector::new_with_value(N, 1.0);

        let a = FVector::new_with_value(N, 2.0);
        let mut b = FVector::new_with_value(N, 3.0);
        b -= &a;
        assert_eq!(b, c);

        let mut b = FVector::new_with_value(N, 3.0);
        b -= a;
        assert_eq!(b, c);
    }
}
