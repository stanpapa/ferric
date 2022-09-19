use crate::linear_algebra::{matrix::FMatrix, matrix_symmetric::FMatrixSym, vector::FVector};

use lapacke::{dgeev, dspev, Layout};

pub trait Diagonalize {
    type Output;

    fn diagonalize(&self) -> Self::Output;
    fn diagonalize_sym(&self) -> Self::Output;
}

// todo: merge FMatrix FMatrixSym impl (macro?)
impl Diagonalize for FMatrix {
    type Output = (FVector, FMatrix);

    fn diagonalize(&self) -> Self::Output {
        if self.cols != self.rows {
            panic!("[diagonalize] trying to diagonalize a non-square matrix.");
        }

        let info;
        let n = self.rows;
        let ldvl = 1;
        let mut a = self.clone();
        let mut eigenvalues_real = FVector::zero(n);
        let mut eigenvalues_imag = FVector::zero(n);
        // let mut eigenvectors_left = FMatrix::new(n, n);
        let mut eigenvectors_left = vec![];
        let mut eigenvectors_right = FMatrix::zero(n, n);

        // todo: not correct
        unsafe {
            info = dgeev(
                Layout::RowMajor,
                b'N',
                b'V',
                n as i32,
                &mut a,
                n as i32,
                &mut eigenvalues_real,
                &mut eigenvalues_imag,
                &mut eigenvectors_left,
                n as i32,
                &mut eigenvectors_right,
                n as i32,
            );
        }

        if info != 0 {
            panic!("Diagonalization failed with error code: {}", info);
        }
        println!("right:\n{}", eigenvectors_right);

        // doesn't give the proper diagonal matrix
        let diag = eigenvectors_right.transposed() * (self.clone() * eigenvectors_right.clone());
        println!("diag:\n{}", diag);

        (eigenvalues_real, eigenvectors_right)
    }

    /// careful: this assumes a symmetric matrix
    fn diagonalize_sym(&self) -> Self::Output {
        let mat_sym = FMatrixSym::from(self);
        mat_sym.diagonalize()
    }
}

impl Diagonalize for FMatrixSym {
    type Output = (FVector, FMatrix);

    fn diagonalize(&self) -> Self::Output {
        self.diagonalize_sym()
    }

    fn diagonalize_sym(&self) -> Self::Output {
        let n = self.n;
        let info;
        let mut ap = self.clone();
        let mut eigenvalues = FVector::zero(n);
        let mut eigenvectors = FMatrix::zero(n, n);

        unsafe {
            info = dspev(
                Layout::RowMajor,
                b'V',
                b'L',
                n as i32,
                &mut ap,
                &mut eigenvalues,
                &mut eigenvectors,
                n as i32,
            );
        }

        if info != 0 {
            panic!("Diagonalization failed with error code: {}", info);
        }

        (eigenvalues, eigenvectors)
    }
}

#[cfg(test)]
mod tests {
    use crate::linear_algebra::diagonalize::Diagonalize;
    use crate::linear_algebra::matrix::FMatrix;
    use crate::linear_algebra::matrix_symmetric::FMatrixSym;

    #[test]
    fn diagonalize() {
        // non-symmetric matrix

        // symmetric matrix
        let a = FMatrix::new_from_vec(3, 3, &[3.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0]);
        let (eigenvalues, _eigenvectors) = a.diagonalize();

        for (one, another) in eigenvalues.iter().zip(&[2.0, 2.0, 5.0]) {
            // 2 5 2 is
            println!("{}, {}", one, another);
            assert!((one - another).abs() < 1e-14);
        }
    }

    #[test]
    fn diagonalize_sym() {
        let a = FMatrix::new_from_vec(3, 3, &[3.0, 1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0]);
        let (eigenvalues, _eigenvectors) = a.diagonalize_sym();

        for (one, another) in eigenvalues.iter().zip(&[2.0, 2.0, 5.0]) {
            println!("{}, {}", one, another);
            assert!((one - another).abs() < 1e-14);
        }
    }

    #[test]
    fn diagonalize_matrix_sym() {
        let a = FMatrixSym::new_from_vec(3, &[3.0, 1.0, 3.0, 1.0, 1.0, 3.0]);
        let (eigenvalues, _eigenvectors) = a.diagonalize();

        for (one, another) in eigenvalues.iter().zip(&[2.0, 2.0, 5.0]) {
            println!("{}, {}", one, another);
            assert!((one - another).abs() < 1e-14);
        }
    }
}
