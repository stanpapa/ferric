use super::matrix::Matrix;
use crate::math::scalar::Scalar;
use std::ops::Mul;

// ---------------------------------------------------------------------
// BLAS Level 3: Matrix Matrix operations
// ---------------------------------------------------------------------

// mat_sum
// trace
// dot
// mat_p_mat
// mat_x_mat
// Hadamard

fn check_mat_mat<T: Scalar>(s: &str, lhs: &Matrix<T>, rhs: &Matrix<T>) {
    if lhs.rows() != rhs.cols() || lhs.cols() != rhs.rows() {
        panic!("[{}] Incompatible dimensions.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

macro_rules! mul_matrix_mat {
    ($($t:ty),*) => {$(
// Matrix = &Matrix * &Matrix -> C(i,j) = A(i,k) B(k,j)
impl Mul<&Matrix<$t>> for &Matrix<$t> {
    type Output = Matrix<$t>;

    fn mul(self, rhs: &Matrix<$t>) -> Self::Output {
        check_mat_mat("Mul", self, rhs);

        let mut mat = Matrix::<$t>::zero(self.rows(), rhs.cols());

        for i in 0..self.rows() {
            for j in 0..rhs.cols() {
                for k in 0..self.cols() {
                    mat[(i, j)] += self[(i, k)] * rhs[(k, j)];
                }
            }
        }

        mat
    }
}
    )*};
}

mul_matrix_mat!(f32, f64, i8, i16, i32, i64, u8, u16, u32, u64);

#[cfg(test)]
mod tests {
    use crate::math::matrix::IMatrix;

    #[test]
    fn mul_matrix_matrix() {
        // mat and vec are initialised in a clunky way on purpose
        // I don't them to be initialized with `vec!` (yet)

        // ( 1  2 )
        // ( 3  2 )
        // ( 1  2 )
        let mut a = IMatrix::new_with_value(3, 2, 2);
        a[(0, 0)] = 1;
        a[(1, 0)] = 3;
        a[(2, 0)] = 1;

        // ( -4  5  6 )
        // (  6 -5  4 )
        let mut b = IMatrix::new_with_value(2, 3, 6);
        b[(0, 0)] = -4;
        b[(0, 1)] = 5;
        b[(1, 1)] = -5;
        b[(1, 2)] = 4;

        // (  8  -5  14 )
        // (  0   5  26 )
        // (  8  -5  14 )
        let mut c = IMatrix::new_with_value(3, 3, 8);
        c[(0, 1)] = -5;
        c[(0, 2)] = 14;
        c[(1, 0)] = 0;
        c[(1, 1)] = 5;
        c[(1, 2)] = 26;
        c[(2, 1)] = -5;
        c[(2, 2)] = 14;

        assert_eq!(&a * &b, c)
    }
}
