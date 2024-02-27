use crate::linear_algebra::{
    diagonalize::{Diagonalize, DiagonalizeSym},
    matrix::FMatrix,
};

/// Calculate power of matrix, non-symmetric variant
pub trait Power {
    type Output;

    fn powi(&self, n: i32) -> Self::Output;
    fn powf(&self, pow: f64) -> Self::Output;
}

/// Calculate power of matrix, symmetric variant
pub trait PowerSym {
    type Output;

    fn powi_sym(&self, n: i32) -> Self::Output;
    fn powf_sym(&self, pow: f64) -> Self::Output;
}

impl Power for FMatrix {
    type Output = FMatrix;

    fn powi(&self, n: i32) -> Self::Output {
        self.powf(f64::from(n))
    }

    fn powf(&self, n: f64) -> Self::Output {
        // diagonalize matrix -> A -> A L = L Λ
        let (eigenvalues, eigenvectors) = self.diagonalize();

        // construct Λ^n
        let mut lambda_n = FMatrix::zero(eigenvalues.n, eigenvalues.n);
        if n == 0.5 {
            for i in 0..lambda_n.cols {
                lambda_n[(i, i)] = eigenvalues[i].sqrt();
            }
        } else if n == -0.5 {
            for i in 0..lambda_n.cols {
                lambda_n[(i, i)] = 1.0 / eigenvalues[i].sqrt();
            }
        } else {
            for i in 0..lambda_n.cols {
                lambda_n[(i, i)] = eigenvalues[i].powf(n);
            }
        }

        // construct A^n = L Λ^n L^T
        &eigenvectors * (lambda_n * eigenvectors.transposed())
    }
}

impl PowerSym for FMatrix {
    type Output = FMatrix;

    fn powi_sym(&self, n: i32) -> Self::Output {
        self.powf_sym(f64::from(n))
    }

    fn powf_sym(&self, n: f64) -> Self::Output {
        // diagonalize matrix -> A -> A L = L Λ
        let (eigenvalues, eigenvectors) = self.diagonalize_sym();

        // construct Λ^n
        let mut lambda_n = FMatrix::zero(eigenvalues.n, eigenvalues.n);
        if n == 0.5 {
            for i in 0..lambda_n.cols {
                lambda_n[(i, i)] = eigenvalues[i].sqrt();
            }
        } else if n == -0.5 {
            for i in 0..lambda_n.cols {
                lambda_n[(i, i)] = 1.0 / eigenvalues[i].sqrt();
            }
        } else {
            for i in 0..lambda_n.cols {
                lambda_n[(i, i)] = eigenvalues[i].powf(n);
            }
        }

        // construct A^n = L Λ^n L^T
        &eigenvectors * (lambda_n * eigenvectors.transposed())
    }
}

// todo!() tests
