use crate::geometry::{distance, Molecule};

pub fn nuclear_repulsion(mol: &Molecule) -> f64 {
    let n = mol.num_atoms();

    let mut vnn = 0.0;
    for i in 0..n {
        for j in 0..i {
            vnn += mol.atoms()[i].z() as f64 * mol.atoms()[j].z() as f64
                / distance(&mol.atoms()[i], &mol.atoms()[j]);
        }
    }

    vnn
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::_ANG_AU;

    #[test]
    fn nuclear_repulsion() {
        let mut mol = Molecule::default();
        mol.scale_coords(&_ANG_AU);

        assert_eq!(super::nuclear_repulsion(&mol), 9.055003146181436);
    }
}
