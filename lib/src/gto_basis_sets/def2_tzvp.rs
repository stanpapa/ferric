use crate::basis::{Basis, Shell};
use crate::geometry::Atom;

pub fn load_def2_tzvp(atoms: &[Atom]) -> Basis {
    println!("Loading def2-tzvp basis set");
    const MAX_ATOMIC_NUMBER: usize = 8;

    let mut shells: Vec<Vec<Shell>> = vec![Default::default(); MAX_ATOMIC_NUMBER + 1];

    // ----------------------------
    // Element #1, Hydrogen
    // ----------------------------
    shells[1] = Vec::with_capacity(4);
    // s
    shells[1].push(Shell::new(
        0,
        vec![34.0613410000, 5.1235746000, 1.1646626000],
        vec![0.0060251978, 0.0450210940, 0.2018972600],
    ));
    // s
    shells[1].push(Shell::new(0, vec![0.3272304100], vec![1.0000000000]));
    // s
    shells[1].push(Shell::new(0, vec![0.1030724100], vec![1.0000000000]));
    // p
    shells[1].push(Shell::new(1, vec![0.8000000000], vec![1.0000000000]));

    // ----------------------------
    // Element #8, Oxygen
    // ----------------------------
    shells[8] = Vec::with_capacity(11);
    // s
    shells[8].push(Shell::new(
        0,
        vec![
            27032.3826310000,
            4052.3871392000,
            922.3272271000,
            261.2407098900,
            85.3546413510,
            31.0350352450,
        ],
        vec![
            0.0002172630,
            0.0016838662,
            0.0087395616,
            0.0352399688,
            0.1115351912,
            0.2558895396,
        ],
    ));
    // s
    shells[8].push(Shell::new(
        0,
        vec![12.2608607280, 4.9987076005],
        vec![0.3976873090, 0.2462784943],
    ));
    // s
    shells[8].push(Shell::new(0, vec![1.1703108158], vec![1.0000000000]));
    // s
    shells[8].push(Shell::new(0, vec![0.4647474099], vec![1.0000000000]));
    // s
    shells[8].push(Shell::new(0, vec![0.1850453636], vec![1.0000000000]));
    // p
    shells[8].push(Shell::new(
        1,
        vec![63.2749548010, 14.6270493790, 4.4501223456, 1.5275799647],
        vec![0.0060685103, 0.0419125758, 0.1615384109, 0.3570695131],
    ));
    // p
    shells[8].push(Shell::new(1, vec![0.5293511794], vec![0.4479420750]));
    // p
    shells[8].push(Shell::new(1, vec![0.1747842127], vec![0.2444606966]));
    // d
    shells[8].push(Shell::new(2, vec![2.3140000000], vec![1.0000000000]));
    // d
    shells[8].push(Shell::new(2, vec![0.6450000000], vec![1.0000000000]));
    // f
    shells[8].push(Shell::new(3, vec![1.4280000000], vec![1.0000000000]));

    Basis::new(atoms, shells)
}
