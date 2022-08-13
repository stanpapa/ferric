use super::{matrix::BaseMatrix, vector::BaseVector};
use crate::tools::math::scalar::Scalar;
use std::ops::Mul;

fn check_mat_vec<T: Scalar>(s: &str, lhs: &BaseMatrix<T>, rhs: &BaseVector<T>) {
    if lhs.cols() != rhs.size() {
        panic!("[{}] Incompatible dimensions.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

fn check_mat_mat<T: Scalar>(s: &str, lhs: &BaseMatrix<T>, rhs: &BaseMatrix<T>) {
    if lhs.rows() != rhs.cols() || lhs.cols() != rhs.rows() {
        panic!("[{}] Incompatible dimensions.", s);
    }

    if !lhs.is_ok() {
        panic!("[{}] lhs is not OK.", s);
    } else if !rhs.is_ok() {
        panic!("[{}] rhs is not OK.", s);
    }
}

// ---------------------------------------------------------------------
// BLAS Level 2
// ---------------------------------------------------------------------

// BaseVector = &BaseMatrix * &BaseVector -> C(i) = A(i,j) B(j)
impl<T: Scalar + Mul<Output = T>> Mul<&BaseVector<T>> for &BaseMatrix<T> {
    type Output = BaseVector<T>;

    fn mul(self, rhs: &BaseVector<T>) -> Self::Output {
        check_mat_vec("Mul", self, rhs);

        let mut vec = BaseVector::<T>::zeros(self.rows());

        for i in 0..self.rows() {
            for j in 0..self.cols() {
                vec[i] += self[(i, j)] * rhs[j];
            }
        }

        vec
    }
}

// ---------------------------------------------------------------------
// BLAS Level 3
// ---------------------------------------------------------------------

// BaseMatrix = &BaseMatrix * &BaseMatrix -> C(i,j) = A(i,k) B(k,j)
impl<T: Scalar + Mul<Output = T>> Mul<&BaseMatrix<T>> for &BaseMatrix<T> {
    type Output = BaseMatrix<T>;

    fn mul(self, rhs: &BaseMatrix<T>) -> Self::Output {
        check_mat_mat("Mul", self, rhs);

        let mut mat = BaseMatrix::<T>::zeros(self.rows(), rhs.cols());

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

#[cfg(test)]
mod tests {
    use crate::tools::math::matrix::IMatrix;
    use crate::tools::math::vector::IVector;

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
