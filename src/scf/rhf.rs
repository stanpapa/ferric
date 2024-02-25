use libferric::geometry::molecule::Molecule;
use libferric::gto_integrals::nuclear_repulsion::nuclear_repulsion;
use libferric::linear_algebra::constants::_AU_EV;
use libferric::linear_algebra::diagonalize::Diagonalize;
use libferric::linear_algebra::matrix_symmetric::FMatrixSym;
use libferric::linear_algebra::power::Power;
use libferric::linear_algebra::traits::Dot;
use libferric::linear_algebra::vector::FVector;
use libferric::linear_algebra::{matrix::FMatrix, matrix_container::FMatrixSymContainer};

use crate::diis::DIIS;
use crate::fock::fock;
use crate::solver::HFSolver;

pub struct RHFSolver {
    c: FMatrix,
    f: FMatrix,
    d: FMatrix,
    eps: FVector,

    homo: usize,
    e: f64,
    nuclear_repulsion: f64,
}

impl RHFSolver {
    pub fn new(c: &FMatrix, molecule: &Molecule) -> Self {
        let homo = molecule.n_electrons / 2;
        println!("HOMO: {:?}", homo);
        let nuclear_repulsion = nuclear_repulsion(molecule.atoms());

        Self {
            c: c.clone(),
            f: FMatrix::new(c.rows, c.cols),
            d: FMatrix::new(c.rows, c.cols),
            eps: FVector::new(c.rows),
            homo,
            e: 0.0,
            nuclear_repulsion,
        }
    }

    fn d_rms(&self, d_old: &FMatrix) -> f64 {
        let mut rms = 0.0;
        let dim = self.d.rows;
        for μ in 0..dim {
            for ν in 0..dim {
                rms += (self.d[(μ, ν)] - d_old[(μ, ν)]).powi(2);
            }
        }
        rms.sqrt()
    }
}

impl HFSolver for RHFSolver {
    fn guess(&mut self, h: &FMatrix, s12: &FMatrix) {
        println!("Guess");
        // --------------------------------
        // build initial guess density
        // --------------------------------
        // core fock matrix: F' = S^{-1/2}^T H S^{-1/2} = S^{-1/2} H S^{-1/2}
        // S-1/2 = symmetric
        // todo implement FMatrixSym * FMatrix
        self.f = s12 * (h * s12);

        // diagonalize F': C'^T F' C' = eps
        let (eps, mut cprime) = self.f.diagonalize_sym();
        self.eps = eps.clone();
        self.c = s12 * cprime;

        // --------------------------------
        // Initial density
        // --------------------------------
        self.density();
    }

    fn density(&mut self) {
        let c_occ = self.c.slice(0, self.c.rows - 1, 0, self.homo - 1);
        self.d = 2.0 * &c_occ * c_occ.transposed();
    }

    fn fock(&mut self, h: &FMatrix, eri: &FMatrixSymContainer) {
        self.f = fock(&self.d, h, eri, 1.0, 0.5);
    }

    fn energy(&mut self, h: &FMatrix) {
        let x = h + self.f.clone();
        self.e = 0.5 * self.d.dot(&x) + self.nuclear_repulsion;
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

        let max_iter = 40;
        let e_thresh = 1e-6;
        let rms_thresh = 1e-12;
        let mut ΔE;
        let mut converged = false;

        let mut diis = DIIS::new(6, 1, s, &s12);
        println!(
            "Iter {:^16} {:^16} {:^16}",
            "E",
            "ΔE",
            "D(rms)" // "Iter {:^16} {:^16} {:^16} {:^16}",
                     // "E", "ΔE", "[D,F]", "D(rms)"
        );
        for iter in 0..max_iter {
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
            if iter >= diis.iter_start {
                diis.do_diis(&mut self.f, &self.d);
            }
            let f_prime = &s12 * &self.f * &s12;
            let (eps, cprime) = f_prime.diagonalize_sym();
            self.eps = eps.clone();
            self.c = &s12 * cprime;
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
            if ΔE.abs() < e_thresh && rms < rms_thresh {
                converged = true;
                break;
            }
        }

        if converged {
            println!("Converged!\n");
        } else {
            println!(
                "Warning: Wavefunction not converged within {} iterations",
                max_iter
            );
        }

        // spin contamination
        self.print_energy(h);
    }

    fn print_energy(&self, h: &FMatrix) {
        println!("----------------");
        println!("Total SCF Energy");
        println!("----------------\n");
        let e1 = h.dot(&self.d);
        let e2 = self.e - e1 - self.nuclear_repulsion;
        println!("                          {:^20}  {:^20}", "Hartree", "eV");
        println!(
            "Total Energy:        {:20.9}  {:20.5}\n",
            self.e,
            self.e * _AU_EV
        );
        println!("Components:");
        println!(
            "Nuclear Repulsion:   {:20.9}  {:20.5}",
            self.nuclear_repulsion,
            self.nuclear_repulsion * _AU_EV
        );
        println!(
            "Electronic Energy:   {:20.9}  {:20.5}",
            self.e - self.nuclear_repulsion,
            (self.e - self.nuclear_repulsion) * _AU_EV
        );
        println!("One Electron Energy: {:20.9}  {:20.5}", e1, e1 * _AU_EV);
        println!("Two Electron Energy: {:20.9}  {:20.5}", e2, e2 * _AU_EV);

        println!("\n\n----------------");
        println!("Orbital Energies");
        println!("----------------\n");
        println!("{}", self.eps);
    }
}
