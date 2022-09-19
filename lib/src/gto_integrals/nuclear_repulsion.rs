use crate::geometry::atom::Atom;
use crate::linear_algebra::functions::distance;

pub fn nuclear_repulsion(mol: &[Atom]) -> f64 {
    let n = mol.len();

    let mut vnn = 0.0;
    for i in 0..n {
        for j in 0..i {
            vnn +=
                mol[i].z() as f64 * mol[j].z() as f64 / distance(mol[i].origin(), mol[j].origin());
        }
    }

    vnn
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::molecule::Molecule;
    use crate::misc::elements::Element::{H, O};

    #[test]
    fn nuclear_repulsion() {
        let mol = Molecule::new(
            vec![
                Atom::new(O, [0.0000000000, 0.0000000000, -0.1190150726]),
                Atom::new(H, [0.7685504811, 0.0000000000, 0.4760602904]),
                Atom::new(H, [-0.7685504811, 0.0000000000, 0.4760602904]),
            ],
            0,
            1,
        );

        assert_eq!(super::nuclear_repulsion(mol.atoms()), 9.055003146181436);
    }
}
