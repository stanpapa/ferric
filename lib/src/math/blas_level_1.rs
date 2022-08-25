// ---------------------------------------------------------------------
// BLAS Level 1: Vector operations
// ---------------------------------------------------------------------

use crate::math::traits::{Dot, Norm};
use crate::math::vector::FVector;
use std::ops::{AddAssign, MulAssign};

use crate::math::matrix::FMatrix;
use crate::math::matrix_symmetric::FMatrixSym;
use crate::math::utils::check_vec_vec;
use blas::*;

// todo:
//   [ ] DSWAP : swap x and y
//   [x] DSCAL : x = a * x
//   [ ] DCOPY : copy x into y
//   [x] DAXPY : y = a * x + y
//   [x] DDOT  : dot product
//   [x] DNRM2 : Euclidean norm

/// Vector Sum
///
// macro_rules! blas_level1 {
//     ($($t:ty),*) => {$(
// impl Vector<$t> {
//
//     fn sum(&self) -> $t {
//         self.iter().sum()
//     }
// }
//
// ///
// /// Vector + Vector
// ///
//
// // ---------------------------------------------------------------------
// // Add
// // ---------------------------------------------------------------------
//
// // Vector = Vector + Vector
// impl Add<Vector<$t>> for Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn add(self, rhs: Vector<$t>) -> Self::Output {
//         Self::add(self, &rhs)
//     }
// }
//
// // Vector = Vector + &Vector
// impl Add<&Vector<$t>> for Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn add(self, rhs: &Vector<$t>) -> Self::Output {
//         check("Add", &self, rhs);
//
//         let mut vec = self.clone();
//         vec += rhs;
//
//         vec
//     }
// }
//
// // Vector = &Vector + Vector
// impl Add<Vector<$t>> for &Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn add(self, rhs: Vector<$t>) -> Self::Output {
//         Self::add(self, &rhs)
//     }
// }
//
// // Vector = &Vector + &Vector
// impl Add<&Vector<$t>> for &Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn add(self, rhs: &Vector<$t>) -> Self::Output {
//         check("Add", self, rhs);
//
//         let mut vec = self.clone();
//         vec += rhs;
//
//         vec
//     }
// }
//
// // ---------------------------------------------------------------------
// // AddAssign
// // ---------------------------------------------------------------------
//

//
// ///
// /// Vector x scalar
// ///
//
// // ---------------------------------------------------------------------
// // Mul
// // ---------------------------------------------------------------------
//
// // Vector = Vector * scalar
// impl Mul<$t> for Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: $t) -> Self::Output {
//         Self::mul(self, &rhs)
//     }
// }
//
// // Vector = Vector * &scalar
// impl Mul<&$t> for Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: &$t) -> Self::Output {
//         let v: Vec<$t> = self.iter().map(|v| v * rhs).collect();
//         Vector::<$t>::new_from_vec(&v[..])
//     }
// }
//
// // Vector = &Vector * scalar
// impl Mul<$t> for &Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: $t) -> Self::Output {
//         Self::mul(self, &rhs)
//     }
// }
//
// // Vector = &Vector * &scalar
// impl Mul<&$t> for &Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: &$t) -> Self::Output {
//         let v: Vec<$t> = self.iter().map(|v| v * rhs).collect();
//         Vector::<$t>::new_from_vec(&v[..])
//     }
// }
//
// // Vector = scalar * Vector
// impl Mul<Vector<$t>> for $t {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: Vector<$t>) -> Self::Output {
//         Self::mul(self, &rhs)
//     }
// }
//
// // Vector = scalar * &Vector
// impl Mul<&Vector<$t>> for $t {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: &Vector<$t>) -> Self::Output {
//         let v: Vec<$t> = rhs.iter().map(|v| v * self).collect();
//         Vector::<$t>::new_from_vec(&v[..])
//     }
// }
//
// // Vector = &scalar * Vector
// impl Mul<Vector<$t>> for &$t {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: Vector<$t>) -> Self::Output {
//         Self::mul(self, &rhs)
//     }
// }
//
// // Vector = &scalar * &Vector
// impl Mul<&Vector<$t>> for &$t {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: &Vector<$t>) -> Self::Output {
//         let v: Vec<$t> = rhs.iter().map(|v| v * self).collect();
//         Vector::<$t>::new_from_vec(&v[..])
//     }
// }
//
// // ---------------------------------------------------------------------
// // MulAssign
// // ---------------------------------------------------------------------
//
// // Vector *= scalar
// impl MulAssign<$t> for Vector<$t> {
//     fn mul_assign(&mut self, rhs: $t) {
//         Self::mul_assign(self, &rhs)
//     }
// }
//
// // Vector *= &scalar
// impl MulAssign<&$t> for Vector<$t> {
//     fn mul_assign(&mut self, rhs: &$t) {
//         for x in self.iter_mut() {
//             *x *= rhs;
//         }
//     }
// }
//
// ///
// /// Vector - Vector
// ///
//
// // ---------------------------------------------------------------------
// // Sub
// // ---------------------------------------------------------------------
//
// // BaseVector<$t> = BaseVector<$t> - BaseVector<$t>
// impl Sub<Vector<$t>> for Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn sub(self, rhs: Vector<$t>) -> Self::Output {
//         Self::sub(self, &rhs)
//     }
// }
//
// // BaseVector<$t> = BaseVector<$t> - &BaseVector<$t>
// impl Sub<&Vector<$t>> for Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn sub(self, rhs: &Vector<$t>) -> Self::Output {
//         check("Sub", &self, rhs);
//
//         let mut vec = self.clone();
//         vec -= rhs;
//
//         vec
//     }
// }
//
// // BaseVector<$t> = &BaseVector<$t> - BaseVector<$t>
// impl Sub<Vector<$t>> for &Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn sub(self, rhs: Vector<$t>) -> Self::Output {
//         Self::sub(self, &rhs)
//     }
// }
//
// // BaseVector<$t> = &BaseVector<$t> - &BaseVector<$t>
// impl Sub<&Vector<$t>> for &Vector<$t> {
//     type Output = Vector<$t>;
//
//     fn sub(self, rhs: &Vector<$t>) -> Self::Output {
//         check("Sub", self, rhs);
//
//         let mut vec = self.clone();
//         vec -= rhs;
//
//         vec
//     }
// }
//
// // ---------------------------------------------------------------------
// // SubAssign
// // ---------------------------------------------------------------------
//
// // allows BaseVector -= BaseVector
// impl SubAssign<Vector<$t>> for Vector<$t> {
//     fn sub_assign(&mut self, rhs: Vector<$t>) {
//         Self::sub_assign(self, &rhs);
//     }
// }
//
// // allows BaseVector -= &BaseVector
// impl SubAssign<&Vector<$t>> for Vector<$t> {
//     fn sub_assign(&mut self, rhs: &Self) {
//         check("Sub", self, rhs);
//
//         for i in 0..self.size() {
//             self[i] -= rhs[i];
//         }
//     }
// }
//
//     ///
//     /// Vector dot product
//     ///
//     impl Vector<$t> {
//         pub fn dot(&self, rhs: &Vector<$t>) -> $t {
//             check("Dot", self, rhs);
//
//             let mut x = <$t>::default();
//
//             for i in 0..self.size() {
//                 x += self[i] * rhs[i];
//             }
//
//             x
//     }
// }
//    )*};
// }
//
// // blas_level1!(f32, f64, i8, i16, i32, i64, u8, u16, u32, u64);
// blas_level1!(f32, i8, i16, i32, i64, u8, u16, u32, u64);

// allows Vector += Vector
impl AddAssign<FVector> for FVector {
    fn add_assign(&mut self, rhs: Self) {
        Self::add_assign(self, &rhs)
    }
}

// allows Vector += &Vector
impl AddAssign<&FVector> for FVector {
    fn add_assign(&mut self, rhs: &FVector) {
        check_vec_vec("AddAssign", self, rhs);

        unsafe {
            daxpy(
                self.size() as i32,
                1.0,
                rhs.as_slice(),
                1,
                self.as_mut_slice(),
                1,
            );
        }
    }
}

impl Dot<f64> for FVector {
    type Output = f64;

    fn dot(&self, rhs: &FVector) -> Self::Output {
        check_vec_vec("Dot", self, rhs);

        unsafe { ddot(self.size() as i32, self.as_slice(), 1, rhs.as_slice(), 1) }
    }
}

impl Norm for FVector {
    type Output = f64;

    fn norm(&self) -> Self::Output {
        unsafe { dnrm2(self.size() as i32, self.as_slice(), 1) }
    }
}

// Vector *= scalar
impl MulAssign<f64> for FVector {
    fn mul_assign(&mut self, alpha: f64) {
        Self::mul_assign(self, &alpha)
    }
}

// Vector *= &scalar
impl MulAssign<&f64> for FVector {
    fn mul_assign(&mut self, alpha: &f64) {
        unsafe {
            dscal(self.size() as i32, *alpha, self.as_mut_slice(), 1);
        }
    }
}

// Matrix *= scalar
impl MulAssign<f64> for FMatrix {
    fn mul_assign(&mut self, alpha: f64) {
        Self::mul_assign(self, &alpha)
    }
}

// Matrix *= &scalar
impl MulAssign<&f64> for FMatrix {
    fn mul_assign(&mut self, alpha: &f64) {
        unsafe {
            dscal(self.size() as i32, *alpha, self.as_mut_slice(), 1);
        }
    }
}

// MatrixSym *= scalar
impl MulAssign<f64> for FMatrixSym {
    fn mul_assign(&mut self, alpha: f64) {
        Self::mul_assign(self, &alpha)
    }
}

// MatrixSym *= &scalar
impl MulAssign<&f64> for FMatrixSym {
    fn mul_assign(&mut self, alpha: &f64) {
        unsafe {
            dscal(self.num_elements() as i32, *alpha, self.as_mut_slice(), 1);
        }
    }
}

impl FVector {
    pub fn axpy(&mut self, alpha: &f64, x: &FVector) {
        check_vec_vec("axpy", self, x);

        unsafe {
            daxpy(
                self.size() as i32,
                *alpha,
                x.as_slice(),
                1,
                self.as_mut_slice(),
                1,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::math::matrix::FMatrix;
    use crate::math::matrix_symmetric::FMatrixSym;
    use crate::math::traits::{Dot, Norm};
    use crate::math::vector::FVector;

    const N: usize = 5;

    //   [ ] DSWAP : swap x and y
    //   [x] DSCAL : x = a * x
    //   [ ] DCOPY : copy x into y
    //   [x] DAXPY : y = a * x + y
    //   [x] DDOT  : dot product
    //   [x] DNRM2 : Euclidean norm

    // #[test]
    // fn sum() {
    //     let a = IVector::new_from_vec(&vec![2, 3, -1]);
    //
    //     assert_eq!(a.sum(), 4);
    // }
    //
    // #[test]
    // fn add() {
    //     let c = FVector::new_with_value(N, 3.0);
    //
    //     let a = FVector::new_with_value(N, 1.0);
    //     let b = FVector::new_with_value(N, 2.0);
    //     assert_eq!(&a + &b, c);
    //     assert_eq!(a + &b, c);
    //
    //     let a = FVector::new_with_value(N, 1.0);
    //     assert_eq!(&a + b, c);
    //
    //     let a = FVector::new_with_value(N, 1.0);
    //     let b = FVector::new_with_value(N, 2.0);
    //     assert_eq!(a + b, c);
    // }

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

    // #[test]
    // fn mul_scalar() {
    //     let vec = FVector::new_with_value(N, 1.0);
    //     assert_eq!(&2.0 * &vec, FVector::new_with_value(N, 2.0));
    //     assert_eq!(2.0 * &vec, FVector::new_with_value(N, 2.0));
    //     assert_eq!(&vec * &2.0, FVector::new_with_value(N, 2.0));
    //     assert_eq!(&vec * 2.0, FVector::new_with_value(N, 2.0));
    //
    //     let vec = FVector::new_with_value(N, 1.0);
    //     assert_eq!(&2.0 * vec, FVector::new_with_value(N, 2.0));
    //
    //     let vec = FVector::new_with_value(N, 1.0);
    //     assert_eq!(2.0 * vec, FVector::new_with_value(N, 2.0));
    //
    //     let vec = FVector::new_with_value(N, 1.0);
    //     assert_eq!(vec * &2.0, FVector::new_with_value(N, 2.0));
    //
    //     let vec = FVector::new_with_value(N, 1.0);
    //     assert_eq!(vec * 2.0, FVector::new_with_value(N, 2.0));
    // }

    #[test]
    fn mul_assign_scalar() {
        let mut vec = FVector::new_with_value(N, 1.0);

        vec *= 2.0;
        assert_eq!(vec, FVector::new_with_value(N, 2.0));

        vec *= &2.0;
        assert_eq!(vec, FVector::new_with_value(N, 4.0));

        let mut mat = FMatrix::new_with_value(N, N, 1.0);
        mat *= 2.0;
        assert_eq!(mat, FMatrix::new_with_value(N, N, 2.0));

        mat *= &2.0;
        assert_eq!(mat, FMatrix::new_with_value(N, N, 4.0));

        let mut mat_sym = FMatrixSym::new_with_value(N, 1.0);
        mat_sym *= 2.0;
        assert_eq!(mat_sym, FMatrixSym::new_with_value(N, 2.0));

        mat_sym *= &2.0;
        assert_eq!(mat_sym, FMatrixSym::new_with_value(N, 4.0));
    }

    // #[test]
    // fn sub() {
    //     let c = FVector::new_with_value(N, 1.0);
    //
    //     let a = FVector::new_with_value(N, 3.0);
    //     let b = FVector::new_with_value(N, 2.0);
    //     assert_eq!(&a - &b, c);
    //     assert_eq!(a - &b, c);
    //
    //     let a = FVector::new_with_value(N, 3.0);
    //     assert_eq!(&a - b, c);
    //
    //     let a = FVector::new_with_value(N, 3.0);
    //     let b = FVector::new_with_value(N, 2.0);
    //     assert_eq!(a - b, c);
    // }
    //
    // #[test]
    // fn sub_assign() {
    //     let c = FVector::new_with_value(N, 1.0);
    //
    //     let a = FVector::new_with_value(N, 2.0);
    //     let mut b = FVector::new_with_value(N, 3.0);
    //     b -= &a;
    //     assert_eq!(b, c);
    //
    //     let mut b = FVector::new_with_value(N, 3.0);
    //     b -= a;
    //     assert_eq!(b, c);
    // }

    #[test]
    fn axpy() {
        let mut y = FVector::new_with_value(N, 2.0);
        let x = FVector::new_with_value(N, 3.0);
        y.axpy(&3.0, &x);

        assert_eq!(y, FVector::new_with_value(N, 11.0));
    }

    #[test]
    fn dot() {
        // let a = IVector::new_from_vec(&vec![2, 3, -1]);
        // let b = IVector::new_from_vec(&vec![-3, 1, 2]);
        //
        // assert_eq!(a.dot(&b), -5);

        let a = FVector::new_from_vec(&vec![2.0, 3.0, -1.0]);
        let b = FVector::new_from_vec(&vec![-3.0, 1.0, 2.0]);

        assert_eq!(a.dot(&b), -5.0);
    }

    #[test]
    fn norm() {
        let a = FVector::new_from_vec(&vec![2.0, 3.0, -1.0]);

        assert_eq!(a.norm(), 14.0_f64.sqrt());
    }
}
