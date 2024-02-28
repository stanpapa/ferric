use crate::{
    geometry::atom::Atom,
    linear_algebra::{
        constants::PI,
        functions::{BinomialCoefficient, Factorial},
        matrix::FMatrix,
    },
    misc::elements::Element,
};

use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fs::File,
    io::prelude::*,
};

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Shell {
    l: u8,
    exps: Vec<f64>,
    coefs: Vec<f64>,
}

impl Shell {
    pub fn new(l: u8, exps: Vec<f64>, coefs: Vec<f64>) -> Self {
        Self { l, exps, coefs }
    }
}

/// return sperical harmonics basis dimension
pub fn dim(l: &u8) -> usize {
    usize::from(2 * l + 1)
}

/// return cartesian basis dimension
pub fn cdim(l: &u8) -> usize {
    usize::from((l + 1) * (l + 2) / 2)
}

#[derive(Default, Debug, Clone, Serialize, Deserialize)]
pub struct BasisShell {
    origin: [f64; 3],
    shell: Shell,

    cbf: Vec<CartesianBasisFunction>,
}

impl BasisShell {
    pub fn new(origin: [f64; 3], shell: Shell) -> Self {
        let mut basis_shell = Self {
            origin,
            shell,
            cbf: vec![],
        };
        basis_shell.generate_cartesian_shell();

        basis_shell
    }

    /// return sperical harmonicas basis dimension
    pub fn dim(&self) -> usize {
        dim(&self.shell.l)
    }

    /// return cartesian basis dimension
    pub fn cdim(&self) -> usize {
        cdim(&self.shell.l)
    }

    // pub fn shell(&self) -> &Shell {
    //     &self.shell
    // }

    pub fn l(&self) -> &u8 {
        &self.shell.l
    }

    pub fn cbf(&self) -> &[CartesianBasisFunction] {
        &self.cbf
    }

    fn generate_cartesian_shell(&mut self) {
        let gaussians: Vec<[u8; 3]> = gaussian_layout(&self.shell.l);
        self.cbf = Vec::with_capacity(gaussians.len());

        for g in gaussians {
            self.cbf.push(CartesianBasisFunction {
                origin: self.origin,
                ml: g,
                exps: self.shell.exps.to_vec(),
                coefs: self.shell.coefs.to_vec(),
                norm: Vec::with_capacity(self.shell.coefs.len()),
            });
        }

        for c in &mut self.cbf {
            c.normalize();
        }
    }
}

#[derive(Default, Clone, Serialize, Deserialize)]
pub struct Basis {
    shells: Vec<BasisShell>,

    l_max: u8,
    cartesian_to_sperical_trafo: HashMap<u8, FMatrix>,
}

impl Basis {
    pub fn new(atoms: &[Atom], shells: Vec<Vec<Shell>>) -> Self {
        let mut basis_shells: Vec<BasisShell> = Default::default();
        let mut l_max = 0;
        for atom in atoms {
            for shell in &shells[atom.z() as usize] {
                basis_shells.push(BasisShell::new(atom.origin.clone(), shell.clone()));
                if shell.l > l_max {
                    l_max = shell.l;
                }
            }
        }

        // construct cartesian to spherical transformatin matrix
        let mut cartesian_to_sperical_trafo = HashMap::<u8, FMatrix>::new();
        for l in 0..=l_max {
            cartesian_to_sperical_trafo.insert(l, Self::cartesian_spherical_transformation(&l));
        }

        Self {
            shells: basis_shells,
            l_max,
            cartesian_to_sperical_trafo,
        }
    }

    /// return sperical harmonicas basis dimension
    pub fn dim(&self) -> usize {
        self.shells.iter().map(|s| s.dim()).sum()
    }

    /// return sperical harmonicas basis dimension
    pub fn cdim(&self) -> usize {
        self.shells.iter().map(|s| s.cdim()).sum()
    }

    pub fn shells(&self) -> &Vec<BasisShell> {
        &self.shells
    }

    pub fn trafo_matrix(&self, l: &u8) -> &FMatrix {
        match self.cartesian_to_sperical_trafo.get(l) {
            Some(mat) => &mat,
            None => panic!("No transformation matrix found for l = {}", l),
        }
    }

    pub fn offset(&self, sn: usize) -> usize {
        (0..sn).map(|i| self.shells[i].dim()).sum()
    }

    // fn coffset(&self, sn: &usize) -> usize {
    //     let mut offset = 0;
    //     for i in 0..*sn {
    //         offset += self.shells[i].cdim();
    //     }
    //     offset
    // }

    /// Cartesian to spherical harmonics transformation matrix
    /// https://doi.org/10.1002/qua.560540202
    fn cartesian_spherical_transformation(lmax: &u8) -> FMatrix {
        let mut exponent: f64;
        let mut s: f64;

        let l = *lmax as i16;

        let xyz = gaussian_layout(lmax);
        let mut mat = FMatrix::zero(dim(lmax), cdim(lmax));

        for c in 0..cdim(lmax) {
            // cartesian components
            let lx = xyz[c][0] as i16;
            let ly = xyz[c][1] as i16;
            let lz = xyz[c][2] as i16;

            for m in -l..=l {
                let mut j = lx + ly - m.abs();
                if j >= 0 && j % 2 == 0 {
                    j /= 2;
                    let mut s1 = 0.0;

                    for i in 0..=((l - m.abs()) / 2) {
                        let mut s2 = 0.0;
                        for k in 0..=j {
                            if (m < 0 && (m.abs() - lx).abs() % 2 == 1)
                                || (m > 0 && (m.abs() - lx).abs() % 2 == 0)
                            {
                                exponent = f64::from((m.abs() - lx + 2 * k) / 2);
                                s = (-1.0_f64).powf(exponent) * 2.0_f64.sqrt();
                            } else if m == 0 && lx % 2 == 0 {
                                exponent = f64::from(k - lx / 2);
                                s = (-1.0f64).powf(exponent);
                            } else {
                                s = 0.0;
                            }

                            s2 += j.binomial_coefficient(&k) as f64
                                * (m.abs().binomial_coefficient(&(lx - 2 * k))) as f64
                                * s;
                        }

                        s1 += l.binomial_coefficient(&i) as f64
                            * i.binomial_coefficient(&j) as f64
                            * (-1.0f64).powf(f64::from(i))
                            * (2 * l - 2 * i).factorial() as f64
                            / (l - m.abs() - 2 * i).factorial() as f64
                            * s2;
                    }

                    mat[((m + l) as usize, c)] = (((2 * lx).factorial() as f64
                        * (2 * ly).factorial() as f64
                        * (2 * lz).factorial() as f64
                        * l.factorial() as f64
                        * (l - m.abs()).factorial() as f64)
                        / (lx.factorial() as f64
                            * ly.factorial() as f64
                            * lz.factorial() as f64
                            * (2 * l).factorial() as f64
                            * (l + m.abs()).factorial() as f64))
                        .sqrt()
                        * s1
                        / (2.0f64.powf(f64::from(l)) * l.factorial() as f64);
                } else {
                    mat[((m + l) as usize, c)] = 0.0;
                }
            }
        }

        mat
    }

    pub fn store(&self, name: &str) {
        let mut buffer = File::create(name.to_owned() + ".basis").expect("Unable to create file");
        write!(
            buffer,
            "{}",
            serde_json::to_string(self).expect("Unable to serialize Basis")
        )
        .expect("Unable to write to file");
    }

    pub fn retrieve(name: &str) -> Self {
        let mut file =
            File::open(name.to_owned() + ".basis").expect("Unable to open file for reading");
        let mut buffer = String::new();
        file.read_to_string(&mut buffer)
            .expect("Unable to read file");
        serde_json::from_str(&buffer).expect("Unable to deserialize a Basis")
    }

    /// Print shell and contraction layout for each unique element present in the calculation
    /// similar to how ORCA prints it
    pub fn print_layout(&self, atoms: &[Atom]) {
        // list of unique elements
        let elements: BTreeSet<Element> = atoms.iter().map(|atom| atom.el.clone()).collect();

        // find a coordinate for each element to read basis set data from
        let origins = elements
            .iter()
            .map(|el| {
                for atom in atoms {
                    if *el == atom.el {
                        return atom.origin;
                    }
                }
                panic!("Could not find an atom corresponding to {}", el);
            })
            .collect::<Vec<[f64; 3]>>();

        // initialise shell layout maps
        // for each element map the angular momentum to the number of primitives
        let mut primitives: BTreeMap<Element, BTreeMap<u8, Vec<u8>>> = elements
            .iter()
            .map(|el| (el.clone(), BTreeMap::<u8, Vec<u8>>::new()))
            .collect();
        // for each element map the angular momentum to the number of contracted functions
        let mut contracted: BTreeMap<Element, BTreeMap<u8, u8>> = elements
            .iter()
            .map(|el| (el.clone(), BTreeMap::<u8, u8>::new()))
            .collect();

        // construct shell layout maps
        let mut count = 0;
        for el in &elements {
            // construct list of basis shells corresponding to the current element
            let shells: Vec<BasisShell> = self
                .shells
                .iter()
                .filter(|shell| shell.origin == origins[count])
                .map(|shell| shell.clone())
                .collect();

            // fill layout
            for shell in shells {
                contracted
                    .get_mut(&el)
                    .unwrap()
                    .entry(*shell.l())
                    .and_modify(|counter| *counter += 1)
                    .or_insert(1);
                primitives
                    .get_mut(&el)
                    .unwrap()
                    .entry(*shell.l())
                    .and_modify(|pattern| pattern.push(shell.shell.coefs.len() as u8))
                    .or_insert(vec![shell.shell.coefs.len() as u8]);
            }

            count += 1;
        }

        // print banner
        println!(
            r#"
---------------------
BASIS SET INFORMATION
---------------------
"#
        );

        // map numerical value of l to orbital label s,p,d,f, etc.
        let labels: HashMap<u8, &str> = vec![(0, "s"), (1, "p"), (2, "d"), (3, "f")]
            .into_iter()
            .collect();

        // for each unique element present shell layout
        for el in &elements {
            println!(
                "{el} {} contracted to {} pattern {{{}}}",
                primitives
                    .get(el)
                    .unwrap()
                    .iter()
                    .map(|(l, pattern)| format!(
                        "{}{}",
                        pattern.iter().map(|&i| i as usize).sum::<usize>(),
                        labels.get(l).unwrap()
                    ))
                    .collect::<String>(),
                contracted
                    .get(el)
                    .unwrap()
                    .iter()
                    .map(|(l, n)| format!("{n}{}", labels.get(l).unwrap()))
                    .collect::<String>(),
                primitives
                    .get(el)
                    .unwrap()
                    .iter()
                    .map(|(_, pattern)| pattern.iter().map(|p| p.to_string()).collect::<String>())
                    .collect::<Vec<String>>()
                    .join("/")
            );
        }
    }

    /// Print basis set information in ORCA format
    pub fn print_orca(&self, atoms: &[Atom]) {
        // list of unique elements
        let elements: BTreeSet<Element> = atoms.iter().map(|atom| atom.el.clone()).collect();

        // map numerical value of l to orbital label s,p,d,f, etc.
        let labels: HashMap<u8, &str> = vec![(0, "S"), (1, "P"), (2, "D"), (3, "F")]
            .into_iter()
            .collect();

        // print banner
        println!(
            r#"
------------------------
BASIS SET IN ORCA FORMAT
------------------------
"#
        );

        // for each unique element present print basis set information
        for element in &elements {
            println!("# Basis set for element : {element}");
            println!("NewGTO {element}");

            // iterate over atoms until we find one that matches the current element
            for atom in atoms {
                if *element == atom.el {
                    for shell in &self.shells {
                        // dangerous to do float comparison?
                        if atom.origin != shell.origin {
                            continue;
                        }

                        // print shell
                        let s = &shell.shell;
                        println!("{} {}", labels.get(&s.l).unwrap(), s.coefs.len());
                        for i in 0..s.coefs.len() {
                            println!(
                                " {:2}  {:17.10} {:17.10}",
                                i + 1,
                                s.exps[i],
                                shell.cbf[0].coefs[i]
                            );
                        }
                    }

                    break;
                }
            }
            println!(" end;\n")
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CartesianBasisFunction {
    origin: [f64; 3],
    ml: [u8; 3],

    exps: Vec<f64>,
    coefs: Vec<f64>,
    norm: Vec<f64>,
}

/// Getters
impl CartesianBasisFunction {
    pub fn origin(&self) -> &[f64; 3] {
        &self.origin
    }

    pub fn ml(&self) -> &[u8; 3] {
        &self.ml
    }

    pub fn ml_i16(&self) -> [i16; 3] {
        [
            i16::from(self.ml[0]),
            i16::from(self.ml[1]),
            i16::from(self.ml[2]),
        ]
    }

    pub fn exps(&self) -> &[f64] {
        &self.exps
    }

    pub fn coefs(&self) -> &[f64] {
        &self.coefs
    }

    pub fn norm(&self) -> &[f64] {
        &self.norm
    }
}

impl CartesianBasisFunction {
    fn normalize(&mut self) {
        let l = self.ml[0];
        let m = self.ml[1];
        let n = self.ml[2];
        let l_total = f64::from(l + m + n);

        // have to cast to signed integer to prevent overflow when
        // either l, m, or n is 0
        let f_ijk: f64 = ((2 * i16::from(l) - 1).factorial2()
            * (2 * i16::from(m) - 1).factorial2()
            * (2 * i16::from(n) - 1).factorial2()) as f64;

        for e in &self.exps {
            let norm_e =
                (2.0_f64.powf(2.0 * l_total + 1.5) * e.powf(l_total + 1.5) / f_ijk / PI.powf(1.5))
                    .sqrt();
            self.norm.push(norm_e);
        }

        let prefactor = PI.powf(1.5) * f_ijk / 2.0_f64.powf(l_total);
        let mut norm = 0.0;

        for ia in 0..self.exps.len() {
            for ib in 0..self.exps.len() {
                norm += self.norm[ia] * self.norm[ib] * self.coefs[ia] * self.coefs[ib]
                    / (self.exps[ia] + self.exps[ib]).powf(l_total + 1.5);
            }
        }

        norm *= prefactor;
        norm = norm.powf(-0.5);

        for c in &mut self.coefs {
            *c *= norm;
        }
    }
}

fn gaussian_layout(l: &u8) -> Vec<[u8; 3]> {
    let n = (l + 1) * (l + 2) / 2;
    let mut layout: Vec<[u8; 3]> = Vec::with_capacity(n as usize);

    let mut a = *l;
    let mut b = 0;
    let mut c = 0;

    let mut count = 0;
    loop {
        layout.push([a, b, c]);

        count += 1;
        if count >= n {
            break;
        }

        if c < l - a {
            b -= 1;
            c += 1;
        } else {
            a -= 1;
            b = l - a;
            c = 0;
        }
    }

    layout
}

#[cfg(test)]
mod tests {
    // #[test]
    // fn test_CBF_normalize() {
    //   let mut cbf = CartesianBasisFunction {
    //     origin: [0.0, 0.0, 0.0],
    //     ml: [1,0,0],
    //     exps: vec![4.2, 2.1, 1.0],
    //     coefs: vec![1.0, 1.0, 1.0],
    //     norm: vec![0.0, 0.0, 0.0],
    //   };
    //   cbf.normalize();
    // }

    use crate::gto_basis_sets::sto_3g;

    #[test]
    fn gaussian_layout() {
        let s = super::gaussian_layout(&0);
        assert_eq!(vec![[0, 0, 0]], s);

        let p = super::gaussian_layout(&1);
        assert_eq!(vec![[1, 0, 0], [0, 1, 0], [0, 0, 1]], p);

        let d = super::gaussian_layout(&2);
        assert_eq!(
            vec![
                [2, 0, 0],
                [1, 1, 0],
                [1, 0, 1],
                [0, 2, 0],
                [0, 1, 1],
                [0, 0, 2]
            ],
            d
        );

        let f = super::gaussian_layout(&3);
        assert_eq!(
            vec![
                [3, 0, 0],
                [2, 1, 0],
                [2, 0, 1],
                [1, 2, 0],
                [1, 1, 1],
                [1, 0, 2],
                [0, 3, 0],
                [0, 2, 1],
                [0, 1, 2],
                [0, 0, 3]
            ],
            f
        );
    }
}
