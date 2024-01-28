use crate::linear_algebra::{
    diagonalize::Diagonalize, matrix::FMatrix, matrix_symmetric::FMatrixSym,
};

pub trait Power {
    type Output;

    fn powi(&self, n: i32) -> Self::Output;
    fn powf(&self, pow: f64) -> Self::Output;
}

impl Power for FMatrixSym {
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
