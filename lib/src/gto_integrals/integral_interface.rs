use crate::basis::Basis;
use crate::geometry::Atom;

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
