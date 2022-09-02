use crate::constants::PI;
use crate::geometry::Atom;
use crate::math::functions::{BinomialCoefficient, Factorial};
use crate::math::matrix::FMatrix;
use std::collections::HashMap;

#[derive(Clone, Default)]
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

/// return sperical harmonicas basis dimension
pub fn dim(l: &u8) -> usize {
    usize::from(2 * l + 1)
}

/// return cartesian basis dimension
pub fn cdim(l: &u8) -> usize {
    usize::from((l + 1) * (l + 2) / 2)
}

#[derive(Default, Clone)]
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

#[derive(Default, Clone)]
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
                basis_shells.push(BasisShell::new(*atom.origin(), shell.clone()));
                if shell.l > l_max {
                    l_max = shell.l;
                }
            }
        }

        // construct cartesian to spherical transformatin matrix
        let mut cartesian_to_sperical_trafo = HashMap::<u8, FMatrix>::new();
        for l in 0..=l_max {
            cartesian_to_sperical_trafo.insert(l, Self::cartesian_spherical_transformation(&l));
            println!("{}", cartesian_to_sperical_trafo.get(&l).unwrap());
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
        let mut offset = 0;
        for i in 0..sn {
            offset += self.shells[i].dim();
        }
        offset
    }

    // fn coffset(&self, sn: &usize) -> usize {
    //     let mut offset = 0;
    //     for i in 0..*sn {
    //         offset += self.shells[i].cdim();
    //     }
    //     offset
    // }

    /// Cartesian to spherical harmenics transformation matrix
    fn cartesian_spherical_transformation(lmax: &u8) -> FMatrix {
        // cartesian components
        let mut lx: u8;
        let mut ly: u8;
        let mut lz: u8;

        let mut exponent: f64;

        let mut s: f64;
        let mut s1: f64;
        let mut s2: f64;

        let l = i16::from(*lmax);

        let xyz = gaussian_layout(lmax);
        let mut mat = FMatrix::zero(dim(lmax), cdim(lmax));

        for c in 0..cdim(lmax) {
            lx = xyz[c][0];
            ly = xyz[c][1];
            lz = xyz[c][2];

            for m in -l..=l {
                let mut j = i16::from(lx) + i16::from(ly) - m.abs();
                if j >= 0 && j % 2 == 0 {
                    j /= 2;
                    s1 = 0.0;

                    for i in 0..=((l - m.abs()) / 2) {
                        s2 = 0.0;
                        for k in 0..=j {
                            if (m < 0 && (m.abs() - i16::from(lx)).abs() % 2 == 1)
                                || (m > 0 && (m.abs() - i16::from(lx)).abs() % 2 == 0)
                            {
                                exponent = f64::from((m.abs() - i16::from(lx) + 2 * k) / 2);
                                s = -1.0_f64.powf(exponent) * 2.0_f64.sqrt();
                            } else if m == 0 && lx % 2 == 0 {
                                exponent = f64::from(k - i16::from(lx) / 2);
                                s = -1.0f64.powf(exponent);
                            } else {
                                s = 0.0;
                            }

                            s2 += j.binomial_coefficient(&k) as f64
                                * (m.abs().binomial_coefficient(&(i16::from(lx) - 2 * k))) as f64
                                * s;
                        }

                        s1 += l.binomial_coefficient(&i) as f64
                            * i.binomial_coefficient(&j) as f64
                            * -1.0f64.powf(f64::from(i))
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
                        / 2.0f64.powf(f64::from(l))
                        * l.factorial() as f64;
                } else {
                    mat[((m + l) as usize, c)] = 0.0;
                }
            }
        }

        mat
    }
}

#[derive(Clone)]
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
