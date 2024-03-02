use crate::geometry::atom::Atom;
use crate::gto_basis_sets::basis::{Basis, Shell};

pub fn load_def2_svp(atoms: &[Atom]) -> Basis {
    println!("Loading def2-SVP basis set");
    const MAX_ATOMIC_NUMBER: usize = 10;

    let mut shells: Vec<Vec<Shell>> = vec![Default::default(); MAX_ATOMIC_NUMBER + 1];


    // ----------------------------
    // Element #1, Hydrogen
    // ----------------------------

    shells[1] = Vec::with_capacity(3);

    // s
    shells[1].push(Shell::new(
        0,
        vec![13.010701, 1.9622572, 0.44453796],
        vec![0.019682158, 0.13796524, 0.47831935],
    ));

    // s
    shells[1].push(Shell::new(
        0,
        vec![0.12194962],
        vec![1.0],
    ));

    // p
    shells[1].push(Shell::new(
        1,
        vec![0.8],
        vec![1.0],
    ));


    // ----------------------------
    // Element #2, Helium
    // ----------------------------

    shells[2] = Vec::with_capacity(3);

    // s
    shells[2].push(Shell::new(
        0,
        vec![38.354936737, 5.7689081479, 1.2399407035],
        vec![0.023814288905, 0.15490906777, 0.46998096633],
    ));

    // s
    shells[2].push(Shell::new(
        0,
        vec![0.29757815953],
        vec![1.0],
    ));

    // p
    shells[2].push(Shell::new(
        1,
        vec![1.0],
        vec![1.0],
    ));


    // ----------------------------
    // Element #3, Lithium
    // ----------------------------

    shells[3] = Vec::with_capacity(5);

    // s
    shells[3].push(Shell::new(
        0,
        vec![266.27785516, 40.069783447, 9.0559944389, 2.4503009051, 0.72209571855],
        vec![0.0064920150325, 0.047747863215, 0.20268796111, 0.48606574817, 0.43626977955],
    ));

    // s
    shells[3].push(Shell::new(
        0,
        vec![0.052810884721],
        vec![1.0],
    ));

    // s
    shells[3].push(Shell::new(
        0,
        vec![0.020960948798],
        vec![1.0],
    ));

    // p
    shells[3].push(Shell::new(
        1,
        vec![1.45, 0.3],
        vec![0.2586, 1.0],
    ));

    // p
    shells[3].push(Shell::new(
        1,
        vec![0.082],
        vec![1.0],
    ));


    // ----------------------------
    // Element #4, Beryllium
    // ----------------------------

    shells[4] = Vec::with_capacity(5);

    // s
    shells[4].push(Shell::new(
        0,
        vec![515.18616125, 77.511037595, 17.552481693, 4.8028940596, 1.4516214316],
        vec![0.0055615307983, 0.041190068062, 0.17913378108, 0.44736716455, 0.4200958192],
    ));

    // s
    shells[4].push(Shell::new(
        0,
        vec![0.1328163362],
        vec![1.0],
    ));

    // s
    shells[4].push(Shell::new(
        0,
        vec![0.045837372213],
        vec![1.0],
    ));

    // p
    shells[4].push(Shell::new(
        1,
        vec![3.6316917145, 0.71695694366, 0.1954193286],
        vec![-0.029033998305, -0.16877854032, -0.51403419628],
    ));

    // p
    shells[4].push(Shell::new(
        1,
        vec![0.06051546589],
        vec![1.0],
    ));


    // ----------------------------
    // Element #5, Boron
    // ----------------------------

    shells[5] = Vec::with_capacity(6);

    // s
    shells[5].push(Shell::new(
        0,
        vec![839.31830086, 126.26464843, 28.620600763, 7.879372271, 2.4088857172],
        vec![-0.0055929201074, -0.041565520771, -0.18299816983, -0.46540391866, -0.44173884791],
    ));

    // s
    shells[5].push(Shell::new(
        0,
        vec![0.25105109036],
        vec![1.0],
    ));

    // s
    shells[5].push(Shell::new(
        0,
        vec![0.083648866069],
        vec![1.0],
    ));

    // p
    shells[5].push(Shell::new(
        1,
        vec![6.0332223619, 1.2499157866, 0.3387167635],
        vec![-0.035603672456, -0.19895775769, -0.50850202618],
    ));

    // p
    shells[5].push(Shell::new(
        1,
        vec![0.096415632351],
        vec![1.0],
    ));

    // d
    shells[5].push(Shell::new(
        2,
        vec![0.5],
        vec![1.0],
    ));


    // ----------------------------
    // Element #6, Carbon
    // ----------------------------

    shells[6] = Vec::with_capacity(6);

    // s
    shells[6].push(Shell::new(
        0,
        vec![1238.4016938, 186.29004992, 42.251176346, 11.676557932, 3.5930506482],
        vec![0.0054568832082, 0.040638409211, 0.18025593888, 0.46315121755, 0.44087173314],
    ));

    // s
    shells[6].push(Shell::new(
        0,
        vec![0.40245147363],
        vec![1.0],
    ));

    // s
    shells[6].push(Shell::new(
        0,
        vec![0.13090182668],
        vec![1.0],
    ));

    // p
    shells[6].push(Shell::new(
        1,
        vec![9.4680970621, 2.0103545142, 0.54771004707],
        vec![0.038387871728, 0.21117025112, 0.51328172114],
    ));

    // p
    shells[6].push(Shell::new(
        1,
        vec![0.15268613795],
        vec![1.0],
    ));

    // d
    shells[6].push(Shell::new(
        2,
        vec![0.8],
        vec![1.0],
    ));


    // ----------------------------
    // Element #7, Nitrogen
    // ----------------------------

    shells[7] = Vec::with_capacity(6);

    // s
    shells[7].push(Shell::new(
        0,
        vec![1712.8415853, 257.64812677, 58.458245853, 16.198367905, 5.0052600809],
        vec![-0.0053934125305, -0.040221581118, -0.1793114499, -0.46376317823, -0.44171422662],
    ));

    // s
    shells[7].push(Shell::new(
        0,
        vec![0.58731856571],
        vec![1.0],
    ));

    // s
    shells[7].push(Shell::new(
        0,
        vec![0.18764592253],
        vec![1.0],
    ));

    // p
    shells[7].push(Shell::new(
        1,
        vec![13.571470233, 2.9257372874, 0.79927750754],
        vec![-0.040072398852, -0.21807045028, -0.51294466049],
    ));

    // p
    shells[7].push(Shell::new(
        1,
        vec![0.21954348034],
        vec![1.0],
    ));

    // d
    shells[7].push(Shell::new(
        2,
        vec![1.0],
        vec![1.0],
    ));


    // ----------------------------
    // Element #8, Oxygen
    // ----------------------------

    shells[8] = Vec::with_capacity(6);

    // s
    shells[8].push(Shell::new(
        0,
        vec![2266.1767785, 340.87010191, 77.363135167, 21.47964494, 6.6589433124],
        vec![-0.0053431809926, -0.03989003923, -0.17853911985, -0.46427684959, -0.44309745172],
    ));

    // s
    shells[8].push(Shell::new(
        0,
        vec![0.80975975668],
        vec![1.0],
    ));

    // s
    shells[8].push(Shell::new(
        0,
        vec![0.25530772234],
        vec![1.0],
    ));

    // p
    shells[8].push(Shell::new(
        1,
        vec![17.721504317, 3.863550544, 1.0480920883],
        vec![0.043394573193, 0.23094120765, 0.51375311064],
    ));

    // p
    shells[8].push(Shell::new(
        1,
        vec![0.27641544411],
        vec![1.0],
    ));

    // d
    shells[8].push(Shell::new(
        2,
        vec![1.2],
        vec![1.0],
    ));


    // ----------------------------
    // Element #9, Fluorine
    // ----------------------------

    shells[9] = Vec::with_capacity(6);

    // s
    shells[9].push(Shell::new(
        0,
        vec![2894.832599, 435.4193912, 98.843328866, 27.485198001, 8.5405498171],
        vec![-0.0053408255515, -0.039904258866, -0.17912768038, -0.46758090825, -0.4465313102],
    ));

    // s
    shells[9].push(Shell::new(
        0,
        vec![1.0654578038],
        vec![1.0],
    ));

    // s
    shells[9].push(Shell::new(
        0,
        vec![0.33247346748],
        vec![1.0],
    ));

    // p
    shells[9].push(Shell::new(
        1,
        vec![22.696633924, 4.9872339257, 1.3491613954],
        vec![-0.045212874436, -0.23754317067, -0.51287353587],
    ));

    // p
    shells[9].push(Shell::new(
        1,
        vec![0.34829881977],
        vec![1.0],
    ));

    // d
    shells[9].push(Shell::new(
        2,
        vec![1.4],
        vec![1.0],
    ));


    // ----------------------------
    // Element #10, Neon
    // ----------------------------

    shells[10] = Vec::with_capacity(6);

    // s
    shells[10].push(Shell::new(
        0,
        vec![3598.9736625, 541.32073112, 122.90450062, 34.216617022, 10.650584124],
        vec![-0.0053259297003, -0.039817417969, -0.17914358188, -0.46893582977, -0.44782537577],
    ));

    // s
    shells[10].push(Shell::new(
        0,
        vec![1.354595396],
        vec![1.0],
    ));

    // s
    shells[10].push(Shell::new(
        0,
        vec![0.41919362639],
        vec![1.0],
    ));

    // p
    shells[10].push(Shell::new(
        1,
        vec![28.424053785, 6.2822510953, 1.6978715079],
        vec![-0.046031944795, -0.23993183041, -0.50871724964],
    ));

    // p
    shells[10].push(Shell::new(
        1,
        vec![0.43300700172],
        vec![1.0],
    ));

    // d
    shells[10].push(Shell::new(
        2,
        vec![1.888],
        vec![1.0],
    ));

    Basis::new(atoms, shells)
}
