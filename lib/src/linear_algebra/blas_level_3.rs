use crate::linear_algebra::{matrix::FMatrix, utils::check_mat_mat};
use std::ops::Mul;

use cblas::{dgemm, Layout, Transpose};

// ---------------------------------------------------------------------
// BLAS Level 3: Matrix Matrix operations
// ---------------------------------------------------------------------

/// todo:
///     [ ] DGEMM - matrix matrix multiply
///     [ ] DSYMM - symmetric matrix matrix multiply

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

        let transa = if transpose_a {
            Transpose::Ordinary
        } else {
            Transpose::None
        };
        let transb = if transpose_b {
            Transpose::Ordinary
        } else {
            Transpose::None
        };

        let m = self.rows as i32;
        let n = self.cols as i32;
        let k = a.cols as i32;
        let lda = if transpose_a { m } else { k };
        let ldb = if transpose_b { k } else { n };

        unsafe {
            dgemm(
                Layout::RowMajor,
                transa,
                transb,
                m,
                n,
                k,
                alpha,
                a,
                lda,
                b,
                ldb,
                beta,
                self,
                n,
            );
        }
    }
}

impl Mul<&FMatrix> for &FMatrix {
    type Output = FMatrix;

    fn mul(self, rhs: &FMatrix) -> Self::Output {
        let mut output = FMatrix::zero(self.rows, rhs.cols);
        output.dgemm(false, false, 1.0, self, rhs, 0.0);
        output
    }
}

impl Mul<&FMatrix> for FMatrix {
    type Output = FMatrix;

    fn mul(self, rhs: &FMatrix) -> Self::Output {
        let mut output = FMatrix::zero(self.rows, rhs.cols);
        output.dgemm(false, false, 1.0, &self, rhs, 0.0);
        output
    }
}

impl Mul<FMatrix> for &FMatrix {
    type Output = FMatrix;

    fn mul(self, rhs: FMatrix) -> Self::Output {
        let mut output = FMatrix::zero(self.rows, rhs.cols);
        output.dgemm(false, false, 1.0, self, &rhs, 0.0);
        output
    }
}

impl Mul<FMatrix> for FMatrix {
    type Output = FMatrix;

    fn mul(self, rhs: FMatrix) -> Self::Output {
        let mut output = FMatrix::zero(self.rows, rhs.cols);
        output.dgemm(false, false, 1.0, &self, &rhs, 0.0);
        output
    }
}

#[cfg(test)]
mod tests {
    use crate::linear_algebra::matrix::FMatrix;

    // #[test]
    // fn dgemm_square_matrix() {}

    #[test]
    fn dgemm_non_square_matrix() {
        let (m, n, k) = (2, 4, 3);
        let a = FMatrix::new_from_vec(m, k, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let b = FMatrix::new_from_vec(
            k,
            n,
            &[
                1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
            ],
        );
        let mut c = FMatrix::new_from_vec(m, n, &[2.0, 6.0, 0.0, 4.0, 7.0, 2.0, 7.0, 2.0]);

        c.dgemm(false, false, 1.0, &a, &b, 1.0);

        assert_eq!(
            c,
            FMatrix::new_from_vec(m, n, &[40.0, 50.0, 50.0, 60.0, 90.0, 100.0, 120.0, 130.0,])
        );
    }

    #[test]
    fn mul() {
        let (m, n, k) = (2, 4, 3);
        let a = FMatrix::new_from_vec(m, k, &[1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        let b = FMatrix::new_from_vec(
            k,
            n,
            &[
                1.0, 5.0, 9.0, 2.0, 6.0, 10.0, 3.0, 7.0, 11.0, 4.0, 8.0, 12.0,
            ],
        );

        assert_eq!(
            a * b,
            FMatrix::new_from_vec(m, n, &[47.0, 53.0, 37.0, 54.0, 89.0, 79.0, 102.0, 103.0,])
        );
    }
}
