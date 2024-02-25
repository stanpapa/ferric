use super::{diis::DIIS, fock::fock, input::SCFInput, solver::HFSolver};

use libferric::{
    geometry::Geometry,
    gto_integrals::nuclear_repulsion::nuclear_repulsion,
    linear_algebra::{
        constants::AU_EV,
        diagonalize::Diagonalize,
        matrix_symmetric::FMatrixSym,
        power::Power,
        traits::Dot,
        vector::FVector,
        {matrix::FMatrix, matrix_container::FMatrixSymContainer},
    },
};

pub struct UHFSolver {
    input: SCFInput,

    c: [FMatrix; 2],
    f: [FMatrix; 2],
    d: [FMatrix; 2],
    eps: [FVector; 2],

    homo: [usize; 2],
    e: f64,
    nuclear_repulsion: f64,
}

impl UHFSolver {
    pub fn new(c: &[FMatrix; 2], geometry: &Geometry, input: SCFInput) -> Self {
        let homo = [geometry.n_electrons_alpha, geometry.n_electrons_beta];
        let nuclear_repulsion = nuclear_repulsion(geometry.molecule.atoms());
        Self {
            input,
            c: [c[0].clone(), c[1].clone()],
            f: [
                FMatrix::new(c[0].rows, c[0].cols),
                FMatrix::new(c[0].rows, c[0].cols),
            ],
            d: [
                FMatrix::new(c[0].rows, c[0].cols),
                FMatrix::new(c[0].rows, c[0].cols),
            ],
            eps: [FVector::new(c[0].rows), FVector::new(c[0].rows)],
            homo,
            e: 0.0,
            nuclear_repulsion,
        }
    }

    fn d_rms(&self, d_old: &[FMatrix; 2]) -> f64 {
        (0..2)
            .map(|op| {
                let mut rms = 0.0;
                let dim = self.d[op].rows;
                for μ in 0..dim {
                    for ν in 0..dim {
                        rms += (self.d[op][(μ, ν)] - d_old[op][(μ, ν)]).powi(2);
                    }
                }
                rms.sqrt()
            })
            .sum()
    }
}

impl HFSolver for UHFSolver {
    fn guess(&mut self, h: &FMatrix, s12: &FMatrix) {
        println!("Guess");
        (0..2).for_each(|op| {
            // --------------------------------
            // build initial guess density
            // --------------------------------
            // core fock matrix: F' = S^{-1/2}^T H S^{-1/2} = S^{-1/2} H S^{-1/2}
            // S-1/2 = symmetric
            // todo implement FMatrixSym * FMatrix
            self.f[op] = s12 * (h * s12);

            // diagonalize F': C'^T F' C' = eps
            let (eps, cprime) = self.f[op].diagonalize_sym();
            self.eps[op] = eps.clone();
            self.c[op] = s12 * cprime;
        });
        // --------------------------------
        // Initial density
        // --------------------------------
        self.density();
    }

    fn density(&mut self) {
        (0..2).for_each(|op: usize| {
            let c_occ = self.c[op].slice(0, self.c[op].rows - 1, 0, self.homo[op] - 1);
            self.d[op] = &c_occ * c_occ.transposed();
        });
    }

    fn fock(&mut self, h: &FMatrix, eri: &FMatrixSymContainer) {
        self.f[0] = fock(&self.d[0], h, eri, 2.0, 1.0);
        self.f[1] = fock(&self.d[1], h, eri, 2.0, 1.0);
    }

    fn energy(&mut self, h: &FMatrix) {
        self.e = (0..2)
            .map(|op: usize| {
                let x = h + self.f[op].clone();
                0.5 * self.d[op].dot(&x)
            })
            .sum::<f64>()
            + self.nuclear_repulsion;
    }
    fn solve(&mut self, h: &FMatrix, eri: &FMatrixSymContainer, s: &FMatrixSym) {
        // --------------------------------
        // build orthogonalization matrix
        // --------------------------------
        let s12 = s.powf(-0.5);

        // --------------------------------
        // Guess
        // --------------------------------
        self.guess(h, &s12);

        let mut ΔE;
        let mut converged = false;

        let mut diis = [
            DIIS::new(self.input.diis_dim_max, self.input.diis_iter_start, s, &s12),
            DIIS::new(6, 1, s, &s12),
        ];
        println!(
            "Iter {:^16} {:^16} {:^16}",
            "E",
            "ΔE",
            "D(rms)" // "Iter {:^16} {:^16} {:^16} {:^16}",
                     // "E", "ΔE", "[D,F]", "D(rms)"
        );
        for iter in 0..self.input.max_iter {
            ΔE = -self.e;
            let d_old = self.d.clone();

            // --------------------------------
            // construct new Fock matrix
            // --------------------------------
            self.fock(h, eri);

            // --------------------------------
            // calculate HF energy
            // --------------------------------
            self.energy(&h);

            // --------------------------------
            // check for convergence
            // --------------------------------
            ΔE += self.e;

            // --------------------------------
            // build new density
            // --------------------------------
            // 1. Orthogonalize F' = S-1/2 F S-1/2
            // 2. Diagonalize F' C' -> C' ε
            // 3. Backtransform: C = S-1/2 C'
            // 4. Compute density: Dμν = \sum_i^{n_occ} Cμi Cνi^T
            if iter >= diis[0].iter_start {
                (0..2).for_each(|op| diis[op].do_diis(&mut self.f[op], &self.d[op]));
            }
            (0..2).for_each(|op| {
                let f_prime = &s12 * &self.f[op] * &s12;
                let (eps, cprime) = f_prime.diagonalize_sym();
                self.eps[op] = eps.clone();
                self.c[op] = &s12 * cprime;
            });
            self.density();

            let rms = self.d_rms(&d_old);
            println!(
                "{:3} {:16.9} {:16.5e} {:16.5e}",
                iter,
                self.e,
                ΔE,
                rms // "{:3} {:16.9} {:16.9} {:16.9} {:16.9}",
                    // iter, e0, ΔE, comm_norm, rms
            );
            if ΔE.abs() < self.input.e_threshold && rms < self.input.rms_threshold {
                converged = true;
                break;
            }
        }

        if converged {
            println!("Converged!\n");
        } else {
            println!(
                "Warning: Wavefunction not converged within {} iterations",
                self.input.max_iter
            );
        }

        self.print_energy(h);
    }

    fn print_energy(&self, h: &FMatrix) {
        println!("----------------");
        println!("Total SCF Energy");
        println!("----------------\n");
        let e1 = (0..2).map(|op| h.dot(&self.d[op])).sum::<f64>();
        let e2 = self.e - e1 - self.nuclear_repulsion;
        println!("                          {:^20}  {:^20}", "Hartree", "eV");
        println!(
            "Total Energy:        {:20.9}  {:20.5}\n",
            self.e,
            self.e * AU_EV
        );
        println!("Components:");
        println!(
            "Nuclear Repulsion:   {:20.9}  {:20.5}",
            self.nuclear_repulsion,
            self.nuclear_repulsion * AU_EV
        );
        println!(
            "Electronic Energy:   {:20.9}  {:20.5}",
            self.e - self.nuclear_repulsion,
            (self.e - self.nuclear_repulsion) * AU_EV
        );
        println!("One Electron Energy: {:20.9}  {:20.5}", e1, e1 * AU_EV);
        println!("Two Electron Energy: {:20.9}  {:20.5}", e2, e2 * AU_EV);

        println!("\n\n----------------");
        println!("Orbital Energies (Alpha)");
        println!("----------------\n");
        println!("{}", self.eps[0]);

        println!("\n\n----------------");
        println!("Orbital Energies (Beta)");
        println!("----------------\n");
        println!("{}", self.eps[1]);
    }
}
