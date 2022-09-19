use crate::{
    geometry::atom::Atom,
    gto_basis_sets::basis::Basis,
    gto_integrals::{one_electron::OneElectronKernel, two_electron::TwoElectronKernel},
};

use std::fs;

pub struct IntegralInterface {
    basis: Basis,
    atoms: Vec<Atom>,
}

impl IntegralInterface {
    pub fn new(basis: &Basis, atoms: &[Atom]) -> IntegralInterface {
        Self {
            basis: basis.clone(),
            atoms: atoms.to_vec(),
        }
    }
}

/// Getters
impl IntegralInterface {
    pub fn basis(&self) -> &Basis {
        &self.basis
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }
}

impl IntegralInterface {
    pub fn remove(&self) {
        for kernel in OneElectronKernel::iter() {
            fs::remove_file(kernel.to_filename()).ok();
        }
        for kernel in TwoElectronKernel::iter() {
            fs::remove_file(kernel.to_filename()).ok();
        }
    }
}
