use super::{matrix::Matrix, vector::BaseVector};
use crate::math::scalar::Scalar;
use std::ops::Mul;

// ---------------------------------------------------------------------
// BLAS Level 2: Matrix Vector operations
// ---------------------------------------------------------------------

// mat_x_vec
// vec_x_mat
// vec_outer_vec

// mat_scal

fn check_mat_vec<T: Scalar>(s: &str, lhs: &Matrix<T>, rhs: &BaseVector<T>) {
    if lhs.cols() != rhs.size() {
        panic!("[{}] Incompatible dimensions.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

// BaseVector = &Matrix * &BaseVector -> C(i) = A(i,j) B(j)
macro_rules! mul_matrix_vec {
    ($($t:ty),*) => {$(
impl Mul<&BaseVector<$t>> for &Matrix<$t> {
    type Output = BaseVector<$t>;

    fn mul(self, rhs: &BaseVector<$t>) -> Self::Output {
        check_mat_vec("Mul", self, rhs);

        let mut vec = BaseVector::<$t>::zero(self.rows());

        for i in 0..self.rows() {
            for j in 0..self.cols() {
                vec[i] += self[(i, j)] * rhs[j];
            }
        }

        vec
    }
}
    )*};
}

mul_matrix_vec!(f32, f64, i8, i16, i32, i64, u8, u16, u32, u64);

#[cfg(test)]
mod tests {
    use crate::math::matrix::IMatrix;
    use crate::math::vector::IVector;

    #[test]
    fn mul_matrix_vector() {
        // mat and vec are initialised in a clunky way on purpose
        // I don't them to be initialized with `vec!` (yet)

        // ( 1  2  3 )
        // ( 4  5  6 )
        // ( 7  8  9 )
        let mut a = IMatrix::new_with_value(3, 3, 1);
        a[(0, 1)] = 2;
        a[(0, 2)] = 3;
        a[(1, 0)] = 4;
        a[(1, 1)] = 5;
        a[(1, 2)] = 6;
        a[(2, 0)] = 7;
        a[(2, 1)] = 8;
        a[(2, 2)] = 9;

        // ( 2  1  3 )
        let mut b = IVector::new_with_value(3, 2);
        b[1] = 1;
        b[2] = 3;

        // ( 13 31 49 )
        let mut c = IVector::new_with_value(3, 13);
        c[1] = 31;
        c[2] = 49;

        assert_eq!(&a * &b, c)
    }
}
