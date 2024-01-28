mod diis;
mod fock;

use crate::diis::DIIS;
use crate::fock::fock;

use ferric_lib::{
    geometry::molecule::Molecule,
    gto_basis_sets::{load_basis_set, BasisSet},
    gto_integrals::{
        integral_interface::IntegralInterface, nuclear_repulsion::nuclear_repulsion,
        one_electron::OneElectronKernel, two_electron::TwoElectronKernel,
    },
    linear_algebra::{
        constants::_AU_EV,
        diagonalize::Diagonalize,
        matrix::FMatrix,
        matrix_container::FMatrixContainer,
        matrix_container::FMatrixSymContainer,
        matrix_symmetric::FMatrixSym,
        power::Power,
        traits::{Dot, Norm},
    },
    HFType,
    HFType::{RHF, UHF},
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

    // --------------------------------
    // build orthogonalization matrix
    // --------------------------------
    let s12 = s.powf(-0.5);

    // --------------------------------
    // build initial guess density
    // --------------------------------
    // core fock matrix: F' = S^{-1/2}^T H S^{-1/2} = S^{-1/2} H S^{-1/2}
    // S-1/2 = symmetric
    // todo implement FMatrixSym * FMatrix
    let h = FMatrix::from(h_core);
    let f = &s12 * (&h * &s12);

    // diagonalize F': C'^T F' C' = eps
    let (mut eps, mut cprime) = f.diagonalize_sym();

    // transform eigenvectors into original (non-orthogonal basis)
    let hf = if molecule.n_electrons_alpha == molecule.n_electrons_beta {
        RHF
    } else {
        UHF
    };
    println!("{} calculation", hf);

    let num_op = if hf == UHF { 2 } else { 1 };
    let mut c = match hf {
        RHF => [&s12 * cprime, FMatrix::default()],
        UHF => {
            let tmp = &s12 * cprime;
            [tmp.clone(), tmp]
        }
        _ => panic!("{} not implemented", hf),
    };

    // construct initial density: Dμν = \sum_i^{n_occ} 2 Cμi Cνi^T
    // retrieve number of electrons and calculate num doubly occupied orbitals
    let homo = match hf {
        RHF => [molecule.n_electrons / 2, 0],
        UHF => [molecule.n_electrons_alpha, molecule.n_electrons_beta],
        _ => panic!("{} not implemented", hf),
    };
    println!("HOMO: {:?}", homo);
    let mut d = [FMatrix::default(), FMatrix::default()];
    for op in 0..num_op {
        d[op] = density(&hf, &c[op], &homo[op]);
    }

    let max_iter = 40;
    let e_thresh = 1e-6;
    let rms_thresh = 1e-12;
    let mut ΔE;
    let mut e0 = 0.0;
    let mut converged = false;

    let mut diis = DIIS::new(6, 1, &s, &s12);
    println!(
        "Iter {:^16} {:^16} {:^16}",
        "E",
        "ΔE",
        "D(rms)" // "Iter {:^16} {:^16} {:^16} {:^16}",
                 // "E", "ΔE", "[D,F]", "D(rms)"
    );
    for iter in 0..max_iter {
        ΔE = -e0;
        let d_old = [d[0].clone(), d[1].clone()];

        // --------------------------------
        // construct new Fock matrix
        // --------------------------------
        let mut f = fock(&h, &d, &eri);

        // --------------------------------
        // calculate HF energy
        // --------------------------------
        e0 = energy(num_op, &d, &h, &f) + nuclear_repulsion;

        // --------------------------------
        // check for convergence
        // --------------------------------
        ΔE += e0;
        // let comm_norm = commutator_d_f(&d[0], &f[0], &s);

        // --------------------------------
        // build new density
        // --------------------------------
        // 1. Orthogonalize F' = S-1/2 F S-1/2
        // 2. Diagonalize F' C' -> C' ε
        // 3. Backtransform: C = S-1/2 C'
        // 4. Compute density: Dμν = \sum_i^{n_occ} Cμi Cνi^T
        match hf {
            RHF => {
                if iter >= diis.iter_start {
                    diis.do_diis(&mut f[0], &d[0]);
                }
                let f_prime = &s12 * &f[0] * &s12;
                (eps, cprime) = f_prime.diagonalize_sym();
                c[0] = &s12 * cprime;
                d[0] = density(&hf, &c[0], &homo[0]);
            }
            UHF => {
                // α
                let f_prime = &s12 * &f[0] * &s12;
                (eps, cprime) = f_prime.diagonalize_sym();
                c[0] = &s12 * cprime;
                d[0] = density(&hf, &c[0], &homo[0]);

                // β
                let f_prime = &s12 * &f[1] * &s12;
                let (_eps, cprime) = f_prime.diagonalize_sym();
                c[1] = &s12 * cprime;
                d[1] = density(&hf, &c[1], &homo[1]);
            }
            _ => panic!("{} not implemented", hf),
        };

        let rms = d_rms(&d[0], &d_old[0]);
        println!(
            "{:3} {:16.9} {:16.9} {:16.9}",
            iter,
            e0,
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
    println!("----------------");
    println!("Total SCF Energy");
    println!("----------------\n");
    let e1 = (0..num_op).map(|op| h.dot(&d[op])).sum::<f64>();
    let e2 = e0 - e1 - nuclear_repulsion;
    println!("                          {:^20}  {:^20}", "Hartree", "eV");
    println!("Total Energy:        {:20.9}  {:20.5}\n", e0, e0 * _AU_EV);
    println!("Components:");
    println!(
        "Nuclear Repulsion:   {:20.9}  {:20.5}",
        nuclear_repulsion,
        nuclear_repulsion * _AU_EV
    );
    println!(
        "Electronic Energy:   {:20.9}  {:20.5}",
        e0 - nuclear_repulsion,
        (e0 - nuclear_repulsion) * _AU_EV
    );
    println!("One Electron Energy: {:20.9}  {:20.5}", e1, e1 * _AU_EV);
    println!("Two Electron Energy: {:20.9}  {:20.5}", e2, e2 * _AU_EV);

    println!("\n\n----------------");
    println!("Orbital Energies");
    println!("----------------\n");
    println!("{}", eps);

    // let c_occ = c[0].slice(0, c[0].rows - 1, 0, homo[0] - 1);
    // println!("Dmo:\n{}", c_occ.transposed() * (&d[0] * &c_occ));

    Ok(())
}

/// Build density as Dμν = \sum_i^{n_occ} Cμi Cνi^T
fn density(hf: &HFType, c: &FMatrix, homo: &usize) -> FMatrix {
    let c_occ = c.slice(0, c.rows - 1, 0, homo - 1);
    match hf {
        RHF => 2.0 * &c_occ * c_occ.transposed(),
        UHF => &c_occ * c_occ.transposed(),
        _ => panic!("[density] {} not implemented", hf),
    }
}

/// Calculate RHF energy as: E0 = 0.5 \sum_{μν} Dμν ( Hμν + Fμν )
fn energy(num_op: usize, d: &[FMatrix; 2], h: &FMatrix, f: &[FMatrix; 2]) -> f64 {
    // let fac = if num_op == 2 { 1.0 } else { 0.5 };
    (0..num_op)
        .map(|op| {
            let x = h + &f[op];
            0.5 * d[op].dot(&x)
        })
        .sum()
}

// /// calculate || [D,F] || with
// ///     [D,F]μν = \sum_{ρσ} Sμρ Dρσ Fσν - Fμρ Dρσ Sσν
// fn commutator_d_f(d: &FMatrix, f: &FMatrix, s: &FMatrixSym) -> f64 {
//     let s_mat = FMatrix::from(s);
//     let mut commutator = &s_mat * (d * f);
//     commutator -= f * (d * s_mat);
//     commutator.norm()
// }

fn d_rms(d: &FMatrix, d_old: &FMatrix) -> f64 {
    let mut rms = 0.0;
    let dim = d.rows;
    for μ in 0..dim {
        for ν in 0..dim {
            rms += (d[(μ, ν)] - d_old[(μ, ν)]).powi(2);
        }
    }
    rms.sqrt()
}
