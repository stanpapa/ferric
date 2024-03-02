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

    pub damp_factor: f64,
    error_max: f64,

    active: bool,

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
            damp_factor: 0.7,
            error_max: f64::MAX,
            active: false,
            error: vec![FMatrix::default(); dim_max],
            fock: vec![FMatrix::default(); dim_max],
        }
    }

    pub fn do_diis(&mut self, f: &mut FMatrix, p: &FMatrix, iter: usize) {
        // don't use DIIS if dim = 0
        if self.dim_max == 0 {
            return;
        }

        // damping becomes counter-productive towards convergence
        let damp_error = 0.1;
        if self.active && iter >= self.iter_start && self.error_max < damp_error {
            self.damp_factor = 0.0;
        }

        // calculate error
        let error = self.calc_error_matrix(f, p);

        // determine max error
        self.error_max = *error
            .iter()
            .max_by(|a, b| a.abs().total_cmp(&b.abs()))
            .unwrap();

        // check start of DIIS
        if !self.active && iter == self.iter_start {
            self.active = true;
            println!("{:^60}", "*** Turning on AO-DIIS ***");
        }

        // in any case store the Fock matrix
        if self.dim < self.dim_max {
            self.dim += 1;
        }

        // store Fock and error matrix
        let pos = iter % self.dim_max;
        self.fock[pos] = f.clone();
        self.error[pos] = error;

        // don't extrapolate if DIIS hasn't been activated yet
        if !self.active {
            return;
        }

        // form DIIS matrix
        let dim = self.dim + 1;
        let mut b = FMatrix::new_with_value(dim, dim, -1.0);
        b[(self.dim, self.dim)] = 0.0;
        for i in 0..self.dim {
            for j in 0..=i {
                let bij = self.error[i].dot(&self.error[j]);
                b[(i, j)] = bij;
                if i != j {
                    b[(j, i)] = bij;
                }
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

    // Error = [P,F] = FPS - SPF
    fn calc_error_matrix(&self, f: &FMatrix, p: &FMatrix) -> FMatrix {
        let mut error = f * (p * &self.s);
        error -= &self.s * (p * f);
        // why do I transform this to an orthonormal basis?
        // &self.s12 * error * &self.s12
        error
    }
}
