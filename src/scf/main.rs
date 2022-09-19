use ferric_lib::{
    geometry::molecule::Molecule,
    gto_basis_sets::{load_basis_set, BasisSet},
    gto_integrals::{
        integral_interface::IntegralInterface, nuclear_repulsion::nuclear_repulsion,
        one_electron::OneElectronKernel, two_electron::TwoElectronKernel,
    },
    linear_algebra::{
        diagonalize::Diagonalize,
        matrix::FMatrix,
        matrix_container::FMatrixContainer,
        matrix_container::FMatrixSymContainer,
        matrix_symmetric::FMatrixSym,
        traits::{Dot, Norm},
    },
};

use std::error;

fn main() -> Result<(), Box<dyn error::Error>> {
    println!("SCF Module");

    // hard-coded base name for now
    let molecule = Molecule::retrieve("input");

    // work-around until I figure out how to store the basis set on disk
    let basis = load_basis_set(&BasisSet::sto_3g, molecule.atoms());
    let integrals = IntegralInterface::new(&basis, molecule.atoms());

    // read integrals from disk
    let nuclear_repulsion = nuclear_repulsion(molecule.atoms());
    let h_core = FMatrixSym::retrieve(OneElectronKernel::HCore.to_filename());
    let s = FMatrixSym::retrieve(OneElectronKernel::Overlap.to_filename());
    let eri = integrals.calc_two_electron_integral(TwoElectronKernel::ERI); // todo: verify
    for keys in eri.keys() {
        println!("{:?}", keys);
    }

    // --------------------------------
    // build orthogonalization matrix
    // --------------------------------
    // diagonalize S -> S L = L Λ
    let (eigenvalues, eigenvectors) = s.diagonalize();

    // construct Λ-1/2
    let mut lambda_12 = FMatrix::zero(eigenvalues.size(), eigenvalues.size());
    for i in 0..lambda_12.cols {
        lambda_12[(i, i)] = 1.0 / eigenvalues[i].sqrt();
    }

    // construct S-1/2 = L Λ-1/2 L^T
    let s12 = &eigenvectors * (lambda_12 * eigenvectors.transposed());
    // println!("S^(-1/2):\n{}", s12); // verified aginst EasyIntegrals

    // --------------------------------
    // build initial guess density
    // --------------------------------
    // core fock matrix: F' = S^{-1/2}^T H S^{-1/2} = S^{-1/2} H S^{-1/2}
    // S-1/2 = symmetric
    // todo implement FMatrixSym * FMatrix
    let h = FMatrix::from(h_core);
    let f = &s12 * (&h * &s12);
    // println!("h:\n{}", h); // verified aginst EasyIntegrals
    // println!("F0:\n{}", f); // verified aginst EasyIntegrals

    // diagonalize F': C'^T F' C' = eps
    let (_eps, cprime) = f.diagonalize_sym();

    // transform eigenvectors into original (non-orthogonal basis)
    let mut c = &s12 * cprime;
    // println!("C0:\n{}", c); // verified aginst EasyIntegrals

    // construct initial density: Dμν = \sum_i^{n_occ} 2 Cμi Cνi^T
    let n_occ = 5; // RHF H20
    let mut d = density(&c, &n_occ);
    // println!("d:\n{}", d); // verified aginst EasyIntegrals

    let max_iter = 20;
    let thresh = 1e-6;
    let mut ΔE;
    let mut e0 = 0.0;
    println!("Iter {:^16} {:^16} {:^16}", "E", "ΔE", "[D,F]");
    for iter in 0..max_iter {
        ΔE = -e0;

        // --------------------------------
        // construct new Fock matrix
        // --------------------------------
        let f = fock(&h, &d, &eri);

        // --------------------------------
        // calculate HF energy
        // --------------------------------
        e0 = energy(&d, &h, &f) + nuclear_repulsion;

        // --------------------------------
        // check for convergence
        // --------------------------------
        ΔE += e0;
        let comm_norm = commutator_d_f(&d, &f, &s);
        println!("{:3} {:16.9} {:16.9} {:16.9}", iter, e0, ΔE, comm_norm);
        if ΔE.abs() < thresh {
            println!("Converged!");
            break;
        }

        // --------------------------------
        // build new density
        // --------------------------------
        // 1. Orthogonalize F' = S-1/2 F S-1/2
        let f_prime = &s12 * &f * &s12;
        // 2. Diagonalize F' -> C'
        let (_eps, cprime) = f_prime.diagonalize_sym();
        // 3. Backtransform: C = S-1/2 C'
        c = &s12 * cprime;
        // 4. Compute density: Dμν = \sum_i^{n_occ} Cμi Cνi^T
        d = density(&c, &n_occ);
    }

    // println!("Dmo:\n{}", c.transposed() * (&d * &c));
    // println!("Fmo:\n{}", c.transposed() * (&f * &c));

    Ok(())
}

/// Build RHF density as Dμν = 2 \sum_i^{n_occ} Cμi Cνi^T
fn density(c: &FMatrix, n: &usize) -> FMatrix {
    let c_occ = c.slice(0, c.rows - 1, 0, n - 1);
    2.0 * &c_occ * c_occ.transposed() // verified against EasyIntegrals
}

/// build RHF Fock matrix as
///     Fμν = Hμν + \sum_{ρσ} Dρσ^{i-1} [ (μν|ρσ) - 0.5 (μρ|νσ) ]
// fn fock(h: &FMatrix, d: &FMatrix, eri: &FMatrixSymContainer) -> FMatrix {
fn fock(h: &FMatrix, d: &FMatrix, eri: &FMatrixContainer) -> FMatrix {
    let mut f = h.clone();
    let dim = f.rows;

    for μ in 0..dim {
        for ν in 0..dim {
            // only lower-triangle is stored
            let μ1 = std::cmp::max(μ, ν);
            let ν1 = std::cmp::min(μ, ν);
            f[(μ, ν)] += d.dot(&eri[(μ1, ν1)]);
        }
    }

    for μ in 0..dim {
        for ν in 0..dim {
            for ρ in 0..dim {
                // only lower-triangle is stored
                let μ1 = std::cmp::max(μ, ρ);
                let ρ1 = std::cmp::min(μ, ρ);
                for σ in 0..dim {
                    f[(μ, ν)] -= 0.5 * d[(ρ, σ)] * eri[(μ1, ρ1)][(ν, σ)];
                }
            }
        }
    }

    f
}

/// Calculate RHF energy as: E0 = 0.5 \sum_{μν} Dμν ( Hμν + Fμν )
fn energy(d: &FMatrix, h: &FMatrix, f: &FMatrix) -> f64 {
    0.5 * d.dot(h) + d.dot(f)
}

/// calculate || [D,F] || with
///     [D,F]μν = \sum_{ρσ} Sμρ Dρσ Fσν - Fμρ Dρσ Sσν
fn commutator_d_f(d: &FMatrix, f: &FMatrix, s: &FMatrixSym) -> f64 {
    let s_mat = FMatrix::from(s);
    let mut commutator = &s_mat * (d * f);
    commutator -= f * (d * s_mat);
    commutator.norm()
}
