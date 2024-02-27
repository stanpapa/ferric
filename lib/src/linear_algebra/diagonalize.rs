use crate::linear_algebra::{matrix::FMatrix, vector::FVector};

use lapacke::{dgeev, dsyev, Layout};

pub trait Diagonalize {
    type Output;

    fn diagonalize(&self) -> Self::Output;
}

pub trait DiagonalizeSym {
    type Output;

    fn diagonalize_sym(&self) -> Self::Output;
}

impl Diagonalize for FMatrix {
    type Output = (FVector, FMatrix);

    fn diagonalize(&self) -> Self::Output {
        if self.cols != self.rows {
            panic!("[diagonalize] trying to diagonalize a non-square matrix.");
        }

        let info;
        let n = self.rows;
        let mut a = self.clone();
        let mut eigenvalues_real = FVector::zero(n);
        let mut eigenvalues_imag = FVector::zero(n);
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
                1,
                &mut eigenvectors_right,
                n as i32,
            );
        }

        if info != 0 {
            panic!("Diagonalization failed with error code: {}", info);
        }

        (eigenvalues_real, eigenvectors_right)
    }
}

impl DiagonalizeSym for FMatrix {
    type Output = (FVector, FMatrix);

    fn diagonalize_sym(&self) -> Self::Output {
        if self.cols != self.rows {
            panic!("[diagonalize] trying to diagonalize a non-square matrix.");
        }

        let info;
        let n = self.rows;
        let mut eigenvectors = self.clone();
        let mut eigenvalues = FVector::zero(n);

        // todo: not correct
        unsafe {
            info = dsyev(
                Layout::RowMajor,
                b'V',
                b'U',
                n as i32,
                &mut eigenvectors,
                n as i32,
                &mut eigenvalues,
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
    use crate::linear_algebra::diagonalize::{Diagonalize, DiagonalizeSym};
    use crate::linear_algebra::matrix::FMatrix;

    #[test]
    fn diagonalize() {
        // non-symmetric matrix
        let a = FMatrix::new_from_vec(3, 3, &[1.0, 2.0, -3.0, 4.0, 5.0, 6.0, -7.0, 8.0, 9.0]);
        let (eigenvalues, _eigenvectors) = a.diagonalize();

        for (one, another) in
            eigenvalues
                .iter()
                .zip(&[-4.749430845008213, 5.215160867742253, 14.534269977265952])
        {
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
}
