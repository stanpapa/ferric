use crate::math::matrix::FMatrix;
use crate::math::matrix_symmetric::FMatrixSym;
use crate::math::utils::{check_mat_sym_vec, check_mat_vec, check_vec_vec};
use crate::math::vector::FVector;

use blas::{dgemv, dspmv, dsymv};

// ---------------------------------------------------------------------
// BLAS Level 2: Matrix Vector operations
// ---------------------------------------------------------------------

/// todo:
///     [x] DGEMV - matrix vector multiply
///     [x] DSPMV - symmetric packed matrix vector multiply
///     [x] DSYMV - symmetric matrix vector multiply
///     [ ] DTRMV - triangular matrix vector multiply

// Vector = &Matrix * &Vector -> C(i) = A(i,j) B(j)
// macro_rules! mul_matrix_vec {
//     ($($t:ty),*) => {$(
// impl Mul<&Vector<$t>> for &Matrix<$t> {
//     type Output = Vector<$t>;
//
//     fn mul(self, rhs: &Vector<$t>) -> Self::Output {
//         check_mat_vec("Mul", self, rhs);
//
//         let mut vec = Vector::<$t>::zero(self.rows());
//
//         for i in 0..self.rows() {
//             for j in 0..self.cols() {
//                 vec[i] += self[(i, j)] * rhs[j];
//             }
//         }
//
//         vec
//     }
// }
//     )*};s
// }
//
// mul_matrix_vec!(f32, f64, i8, i16, i32, i64, u8, u16, u32, u64);

impl FVector {
    // y = alpha * A * x + beta * y
    pub fn dgemv(&mut self, transpose: bool, alpha: f64, a: &FMatrix, x: &FVector, beta: f64) {
        check_mat_vec("dgemv", a, x);
        if a.rows() != self.size() {
            panic!("[dgemv] Matrix Vector product not compatible with y");
        }

        let trans = {
            if transpose {
                b'T'
            } else {
                b'N'
            }
        };

        unsafe {
            dgemv(
                trans,
                a.rows() as i32,
                a.cols() as i32,
                alpha,
                a.as_slice(),
                a.rows() as i32,
                x.as_slice(),
                1,
                beta,
                self.as_mut_slice(),
                1,
            );
        }
    }

    // y = alpha * A * x + beta * y
    pub fn dspmv(&mut self, alpha: f64, a: &FMatrixSym, x: &FVector, beta: f64) {
        check_vec_vec("dspmv", x, self);
        check_mat_sym_vec("dspmv", a, x);

        unsafe {
            dspmv(
                b'L',
                a.rows() as i32,
                alpha,
                a.as_slice(),
                x.as_slice(),
                1,
                beta,
                self.as_mut_slice(),
                1,
            );
        }
    }

    /// y = alpha * A * x + beta * y
    /// Inferior to dspmv. Implemented for completeness' sake
    pub fn dsymv(&mut self, alpha: f64, a: &FMatrixSym, x: &FVector, beta: f64) {
        check_vec_vec("dsymv", x, self);
        check_mat_sym_vec("dsymv", a, x);

        let a_sym = FMatrix::from(a);

        unsafe {
            dsymv(
                b'L',
                a_sym.rows() as i32,
                alpha,
                a_sym.as_slice(),
                a.rows() as i32,
                x.as_slice(),
                1,
                beta,
                self.as_mut_slice(),
                1,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    const M: usize = 4;
    const N: usize = 3;

    // #[test]
    // fn mul_matrix_vector() {
    //     // mat and vec are initialised in a clunky way on purpose
    //     // I don't them to be initialized with `vec!` (yet)
    //
    //     // ( 1  2  3 )
    //     // ( 4  5  6 )
    //     // ( 7  8  9 )
    //     let mut a = IMatrix::new_with_value(3, 3, 1);
    //     a[(0, 1)] = 2;
    //     a[(0, 2)] = 3;
    //     a[(1, 0)] = 4;
    //     a[(1, 1)] = 5;
    //     a[(1, 2)] = 6;
    //     a[(2, 0)] = 7;
    //     a[(2, 1)] = 8;
    //     a[(2, 2)] = 9;
    //
    //     // ( 2  1  3 )
    //     let mut b = IVector::new_with_value(3, 2);
    //     b[1] = 1;
    //     b[2] = 3;
    //
    //     // ( 13 31 49 )
    //     let mut c = IVector::new_with_value(3, 13);
    //     c[1] = 31;
    //     c[2] = 49;
    //
    //     assert_eq!(&a * &b, c)
    // }

    use crate::math::matrix::FMatrix;
    use crate::math::matrix_symmetric::FMatrixSym;
    use crate::math::vector::FVector;

    #[test]
    fn dgemv() {
        // todo: a should have all different values
        let a = FMatrix::new_with_value(M, N, 2.0);
        let x = FVector::new_from_vec(&[-2.0, 3.0, -4.0]);
        let mut y = FVector::new_from_vec(&[1.0, 2.0, 3.0, 4.0]);
        y.dgemv(false, 2.0, &a, &x, 5.0);

        assert_eq!(y, FVector::new_from_vec(&[-7.0, -2.0, 3.0, 8.0]));
    }

    #[test]
    fn dspmv() {
        // todo: a should have all different values
        // todo: fix: matrix vector product is wrong
        let a = FMatrixSym::new_with_value(N, 2.0);
        let x = FVector::new_from_vec(&[-2.0, 3.0, -4.0]);
        let mut y = FVector::new_from_vec(&[1.0, 2.0, 3.0]);
        y.dspmv(2.0, &a, &x, 5.0);

        assert_eq!(y, FVector::new_from_vec(&[-7.0, -2.0, 3.0]));
    }

    #[test]
    fn dsymv() {
        // todo: a should have all different values
        // todo: fix: matrix vector product is wrong
        let a = FMatrixSym::new_with_value(N, 2.0);
        let x = FVector::new_from_vec(&[-2.0, 3.0, -4.0]);
        let mut y = FVector::new_from_vec(&[1.0, 2.0, 3.0]);
        y.dsymv(2.0, &a, &x, 5.0);

        assert_eq!(y, FVector::new_from_vec(&[-7.0, -2.0, 3.0]));
    }
}
