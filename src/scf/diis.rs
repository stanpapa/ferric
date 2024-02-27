use libferric::linear_algebra::{
    linear_solve::LinearSolve, matrix::FMatrix, traits::Dot, vector::FVector,
};

use std::cmp::Ordering;

pub struct DIIS {
    dim_max: usize,
    dim: usize,
    pub iter_start: usize,
    iter_current: usize,

    s: FMatrix,
    s12: FMatrix,

    error: Vec<FMatrix>,
    fock: Vec<FMatrix>,
}

impl DIIS {
    pub fn new(dim_max: usize, iter_start: usize, s: &FMatrix, s12: &FMatrix) -> Self {
        Self {
            dim_max,
            dim: 0,
            iter_start,
            iter_current: 0,
            s: s.clone(),
            s12: s12.clone(),
            error: vec![FMatrix::default(); dim_max],
            fock: vec![FMatrix::default(); dim_max],
        }
    }

    pub fn do_diis(&mut self, f: &mut FMatrix, d: &FMatrix) {
        self.iter_current += 1;

        match self.iter_current.cmp(&self.iter_start) {
            Ordering::Less => return,
            Ordering::Equal => println!("{:^60}", "*** Turning on DIIS ***"),
            _ => (),
        };

        let pos = (self.iter_current - self.iter_start) % (self.dim_max);

        self.fock[pos] = f.clone();
        self.error[pos] = self.calc_error_matrix(&self.fock[pos], d);

        if self.dim < self.dim_max {
            self.dim += 1;
        }

        if self.iter_current <= self.iter_start {
            return;
        }

        let dim = self.dim + 1;
        let mut b = FMatrix::new_with_value(dim, dim, -1.0);
        b[(self.dim, self.dim)] = 0.0;
        for i in 0..self.dim {
            for j in 0..self.dim {
                b[(i, j)] = self.error[i].dot(&self.error[j]);
            }
        }

        // solve linear equation
        let mut rhs = FVector::zero(dim);
        rhs[self.dim] = -1.0;
        b.linear_solve_mut(&mut rhs);

        // construct new guess: F' = \sum_i c_i F_i
        f.init(0.0);
        for i in 0..self.dim {
            *f += rhs[i] * &self.fock[i];
        }
    }

    // Error = FDS - SDF
    fn calc_error_matrix(&self, f: &FMatrix, d: &FMatrix) -> FMatrix {
        let mut error = f * (d * &self.s);
        error -= &self.s * (d * f);
        &self.s12 * error * &self.s12
        // error
    }
}
