use crate::basis::{Basis, Shell};
use crate::geometry::Atom;

pub fn load_sto_3g(atoms: &[Atom]) -> Basis {
    println!("Loading STO-3G basis set");
    const MAX_ATOMIC_NUMBER: usize = 8;

    let mut shells: Vec<Vec<Shell>> = vec![Default::default(); MAX_ATOMIC_NUMBER + 1];

    // ----------------------------
    // Element #1, Hydrogen
    // ----------------------------
    shells[1] = Vec::with_capacity(1);
    // s
    shells[1].push(Shell::new(
        0,
        vec![3.42525091, 0.62391373, 0.16885540],
        vec![0.15432897, 0.53532814, 0.44463454],
    ));

    // ----------------------------
    // Element #8, Oxygen
    // ----------------------------
    shells[8] = Vec::with_capacity(3);
    // s
    shells[8].push(Shell::new(
        0,
        vec![130.7093200000, 23.8088610000, 6.4436083000],
        vec![0.1543289687, 0.5353281356, 0.4446345363],
    ));
    // s
    shells[8].push(Shell::new(
        0,
        vec![5.0331513000, 1.1695961000, 0.3803890000],
        vec![-0.0999672287, 0.3995128246, 0.7001154606],
    ));
    // p
    shells[8].push(Shell::new(
        1,
        vec![5.0331513000, 1.1695961000, 0.3803890000],
        vec![0.1559162685, 0.6076837141, 0.3919573862],
    ));

    Basis::new(atoms, shells)
}
