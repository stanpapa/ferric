// ---------------------------------------------------------------------
// BLAS Level 1: Vector operations
// ---------------------------------------------------------------------

use crate::math::scalar::Scalar;
use crate::math::vector::Vector;
use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign};

/// todo:
///  - vec + scal
///  - Hadamard
//   - vec_norm
//   - mat * scal

fn check<T: Scalar>(s: &str, lhs: &Vector<T>, rhs: &Vector<T>) {
    if lhs.size() != rhs.size() {
        panic!("[{}] Vector dimensions are not the same.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

///
/// Vector Sum
///
macro_rules! blas_level1 {
    ($($t:ty),*) => {$(
impl Vector<$t> {

    fn sum(&self) -> $t {
        self.iter().sum()
    }
}

///
/// Vector + Vector
///

// ---------------------------------------------------------------------
// Add
// ---------------------------------------------------------------------

// Vector = Vector + Vector
impl Add<Vector<$t>> for Vector<$t> {
    type Output = Vector<$t>;

    fn add(self, rhs: Vector<$t>) -> Self::Output {
        Self::add(self, &rhs)
    }
}

// Vector = Vector + &Vector
impl Add<&Vector<$t>> for Vector<$t> {
    type Output = Vector<$t>;

    fn add(self, rhs: &Vector<$t>) -> Self::Output {
        check("Add", &self, rhs);

        let mut vec = self.clone();
        vec += rhs;

        vec
    }
}

// Vector = &Vector + Vector
impl Add<Vector<$t>> for &Vector<$t> {
    type Output = Vector<$t>;

    fn add(self, rhs: Vector<$t>) -> Self::Output {
        Self::add(self, &rhs)
    }
}

// Vector = &Vector + &Vector
impl Add<&Vector<$t>> for &Vector<$t> {
    type Output = Vector<$t>;

    fn add(self, rhs: &Vector<$t>) -> Self::Output {
        check("Add", self, rhs);

        let mut vec = self.clone();
        vec += rhs;

        vec
    }
}

// ---------------------------------------------------------------------
// AddAssign
// ---------------------------------------------------------------------

// allows Vector += Vector
impl AddAssign<Vector<$t>> for Vector<$t> {
    fn add_assign(&mut self, rhs: Self) {
        Self::add_assign(self, &rhs)
    }
}

// allows Vector += &Vector
impl AddAssign<&Vector<$t>> for Vector<$t> {
    fn add_assign(&mut self, rhs: &Self) {
        check("AddAssign", self, rhs);

        for i in 0..self.size() {
            self[i] += rhs[i];
        }
    }
}

///
/// Vector x scalar
///

// ---------------------------------------------------------------------
// Mul
// ---------------------------------------------------------------------

// Vector = Vector * scalar
impl Mul<$t> for Vector<$t> {
    type Output = Vector<$t>;

    fn mul(self, rhs: $t) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// Vector = Vector * &scalar
impl Mul<&$t> for Vector<$t> {
    type Output = Vector<$t>;

    fn mul(self, rhs: &$t) -> Self::Output {
        let v: Vec<$t> = self.iter().map(|v| v * rhs).collect();
        Vector::<$t>::new_from_vec(&v[..])
    }
}

// Vector = &Vector * scalar
impl Mul<$t> for &Vector<$t> {
    type Output = Vector<$t>;

    fn mul(self, rhs: $t) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// Vector = &Vector * &scalar
impl Mul<&$t> for &Vector<$t> {
    type Output = Vector<$t>;

    fn mul(self, rhs: &$t) -> Self::Output {
        let v: Vec<$t> = self.iter().map(|v| v * rhs).collect();
        Vector::<$t>::new_from_vec(&v[..])
    }
}

// Vector = scalar * Vector
impl Mul<Vector<$t>> for $t {
    type Output = Vector<$t>;

    fn mul(self, rhs: Vector<$t>) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// Vector = scalar * &Vector
impl Mul<&Vector<$t>> for $t {
    type Output = Vector<$t>;

    fn mul(self, rhs: &Vector<$t>) -> Self::Output {
        let v: Vec<$t> = rhs.iter().map(|v| v * self).collect();
        Vector::<$t>::new_from_vec(&v[..])
    }
}

// Vector = &scalar * Vector
impl Mul<Vector<$t>> for &$t {
    type Output = Vector<$t>;

    fn mul(self, rhs: Vector<$t>) -> Self::Output {
        Self::mul(self, &rhs)
    }
}

// Vector = &scalar * &Vector
impl Mul<&Vector<$t>> for &$t {
    type Output = Vector<$t>;

    fn mul(self, rhs: &Vector<$t>) -> Self::Output {
        let v: Vec<$t> = rhs.iter().map(|v| v * self).collect();
        Vector::<$t>::new_from_vec(&v[..])
    }
}

// ---------------------------------------------------------------------
// MulAssign
// ---------------------------------------------------------------------

// Vector *= scalar
impl MulAssign<$t> for Vector<$t> {
    fn mul_assign(&mut self, rhs: $t) {
        Self::mul_assign(self, &rhs)
    }
}

// Vector *= &scalar
impl MulAssign<&$t> for Vector<$t> {
    fn mul_assign(&mut self, rhs: &$t) {
        for x in self.iter_mut() {
            *x *= rhs;
        }
    }
}

///
/// Vector - Vector
///

// ---------------------------------------------------------------------
// Sub
// ---------------------------------------------------------------------

// BaseVector<$t> = BaseVector<$t> - BaseVector<$t>
impl Sub<Vector<$t>> for Vector<$t> {
    type Output = Vector<$t>;

    fn sub(self, rhs: Vector<$t>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// BaseVector<$t> = BaseVector<$t> - &BaseVector<$t>
impl Sub<&Vector<$t>> for Vector<$t> {
    type Output = Vector<$t>;

    fn sub(self, rhs: &Vector<$t>) -> Self::Output {
        check("Sub", &self, rhs);

        let mut vec = self.clone();
        vec -= rhs;

        vec
    }
}

// BaseVector<$t> = &BaseVector<$t> - BaseVector<$t>
impl Sub<Vector<$t>> for &Vector<$t> {
    type Output = Vector<$t>;

    fn sub(self, rhs: Vector<$t>) -> Self::Output {
        Self::sub(self, &rhs)
    }
}

// BaseVector<$t> = &BaseVector<$t> - &BaseVector<$t>
impl Sub<&Vector<$t>> for &Vector<$t> {
    type Output = Vector<$t>;

    fn sub(self, rhs: &Vector<$t>) -> Self::Output {
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
impl SubAssign<Vector<$t>> for Vector<$t> {
    fn sub_assign(&mut self, rhs: Vector<$t>) {
        Self::sub_assign(self, &rhs);
    }
}

// allows BaseVector -= &BaseVector
impl SubAssign<&Vector<$t>> for Vector<$t> {
    fn sub_assign(&mut self, rhs: &Self) {
        check("Sub", self, rhs);

        for i in 0..self.size() {
            self[i] -= rhs[i];
        }
    }
}

    ///
    /// Vector dot product
    ///
    impl Vector<$t> {
        pub fn dot(&self, rhs: &Vector<$t>) -> $t {
            check("Dot", self, rhs);

            let mut x = <$t>::default();

            for i in 0..self.size() {
                x += self[i] * rhs[i];
            }

            x
    }
}
   )*};
}

blas_level1!(f32, f64, i8, i16, i32, i64, u8, u16, u32, u64);

#[cfg(test)]
mod tests {
    use crate::math::vector::{FVector, IVector};

    const N: usize = 5;

    #[test]
    fn sum() {
        let a = IVector::new_from_vec(&vec![2, 3, -1]);

        assert_eq!(a.sum(), 4);
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

    #[test]
    fn dot() {
        let a = IVector::new_from_vec(&vec![2, 3, -1]);
        let b = IVector::new_from_vec(&vec![-3, 1, 2]);

        assert_eq!(a.dot(&b), -5);
    }
}
