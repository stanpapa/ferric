use crate::math::matrix::FMatrix;
use crate::math::utils::check_mat_mat;
use blas::dgemm;

// ---------------------------------------------------------------------
// BLAS Level 3: Matrix Matrix operations
// ---------------------------------------------------------------------

/// todo:
///     [ ] DGEMM - matrix matrix multiply
///     [ ] DSYMM - symmetric matrix matrix multiply

// macro_rules! mul_matrix_mat {
//     ($($t:ty),*) => {$(
// // Matrix = &Matrix * &Matrix -> C(i,j) = A(i,k) B(k,j)
// impl Mul<&Matrix<$t>> for &Matrix<$t> {
//     type Output = Matrix<$t>;
//
//     fn mul(self, rhs: &Matrix<$t>) -> Self::Output {
//         check_mat_mat("Mul", self, rhs);
//
//         let mut mat = Matrix::<$t>::zero(self.rows(), rhs.cols());
//
//         for i in 0..self.rows() {
//             for j in 0..rhs.cols() {
//                 for k in 0..self.cols() {
//                     mat[(i, j)] += self[(i, k)] * rhs[(k, j)];
//                 }
//             }
//         }
//
//         mat
//     }
// }
//     )*};
// }
//
// mul_matrix_mat!(f32, f64, i8, i16, i32, i64, u8, u16, u32, u64);

impl FMatrix {
    // C(m,n) = alpha * A(m,k) * B(k,n) + beta * C(m,n)
    pub fn dgemm(
        &mut self,
        transpose_a: bool,
        transpose_b: bool,
        alpha: f64,
        a: &FMatrix,
        b: &FMatrix,
        beta: f64,
    ) {
        check_mat_mat("dgemm", self, a, b);

        let trans_a = {
            if transpose_a {
                b'T'
            } else {
                b'N'
            }
        };
        let trans_b = {
            if transpose_b {
                b'T'
            } else {
                b'N'
            }
        };

        let m = self.rows() as i32;
        let n = self.cols() as i32;
        let k = a.cols() as i32;

        unsafe {
            dgemm(
                trans_a,
                trans_b,
                m,
                n,
                k,
                alpha,
                a.as_slice(),
                m,
                b.as_slice(),
                k,
                beta,
                self.as_mut_slice(),
                m,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::math::matrix::FMatrix;

    // #[test]
    // fn mul_matrix_matrix() {
    //     // mat and vec are initialised in a clunky way on purpose
    //     // I don't them to be initialized with `vec!` (yet)
    //
    //     // ( 1  2 )
    //     // ( 3  2 )
    //     // ( 1  2 )
    //     let mut a = IMatrix::new_with_value(3, 2, 2);
    //     a[(0, 0)] = 1;
    //     a[(1, 0)] = 3;
    //     a[(2, 0)] = 1;
    //
    //     // ( -4  5  6 )
    //     // (  6 -5  4 )
    //     let mut b = IMatrix::new_with_value(2, 3, 6);
    //     b[(0, 0)] = -4;
    //     b[(0, 1)] = 5;
    //     b[(1, 1)] = -5;
    //     b[(1, 2)] = 4;
    //
    //     // (  8  -5  14 )
    //     // (  0   5  26 )
    //     // (  8  -5  14 )
    //     let mut c = IMatrix::new_with_value(3, 3, 8);
    //     c[(0, 1)] = -5;
    //     c[(0, 2)] = 14;
    //     c[(1, 0)] = 0;
    //     c[(1, 1)] = 5;
    //     c[(1, 2)] = 26;
    //     c[(2, 1)] = -5;
    //     c[(2, 2)] = 14;
    //
    //     assert_eq!(&a * &b, c)
    // }

    #[test]
    fn dgemm() {
        let (m, n, k) = (2, 4, 3);
        let a = FMatrix::new_from_vec(m, k, &[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        let b = FMatrix::new_from_vec(
            k,
            n,
            &[
                1.0, 5.0, 9.0, 2.0, 6.0, 10.0, 3.0, 7.0, 11.0, 4.0, 8.0, 12.0,
            ],
        );
        let mut c = FMatrix::new_from_vec(m, n, &[2.0, 7.0, 6.0, 2.0, 0.0, 7.0, 4.0, 2.0]);

        c.dgemm(false, false, 1.0, &a, &b, 1.0);

        assert!(
            c == FMatrix::new_from_vec(m, n, &[40.0, 90.0, 50.0, 100.0, 50.0, 120.0, 60.0, 130.0,])
        );
    }
}
