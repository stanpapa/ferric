use crate::linear_algebra::{
    matrix::FMatrix,
    vector::{FVector, Vector},
};

use lapacke::{dgesv, dspsv, dsysv, Layout};

/// Solve linear system of equations
///        A x = b
pub trait LinearSolve<T> {
    type Output;

    fn linear_solve(&self, b: &T) -> Self::Output;
    fn linear_solve_mut(&self, b: &mut T);
}

impl LinearSolve<FVector> for FMatrix {
    type Output = FVector;

    fn linear_solve(&self, b: &FVector) -> Self::Output {
        let mut x = b.clone();
        self.linear_solve_mut(&mut x);
        x
    }

    fn linear_solve_mut(&self, b: &mut FVector) {
        // check dimensions
        if self.rows != self.cols {
            panic!(
                "[dgesv]: A needs to be square! Current shape is A[{},{}]",
                self.rows, self.cols
            );
        }
        if self.cols != b.n {
            panic!(
                "[dgesv] A and B dimension do not match. A[{},{}], B[{}]",
                self.rows, self.cols, b.n
            );
        }

        let n = b.n;
        let mut tmp = self.clone();

        let mut ipiv = Vector::<i32>::zero(n); // ????
        let info;

        // not working correctly! Numpy works fine
        unsafe {
            info = dgesv(
                Layout::RowMajor,
                n as i32,
                1,
                &mut tmp,
                n as i32,
                &mut ipiv,
                b,
                1,
            );
        }

        if info != 0 {
            panic!("[dsysv] Failed to solve Ax=B, error code {}", info);
        }
    }
}
