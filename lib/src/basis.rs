// use super::gto_bases::*;
use crate::constants::PI;
use crate::geometry::Atom;
use crate::math::functions::Factorial;
use crate::math::matrix::{FMatrix, Matrix};

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
pub fn sdim(l: &u8) -> usize {
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
    pub fn sdim(&self) -> usize {
        sdim(&self.shell.l)
    }

    /// return cartesian basis dimension
    pub fn cdim(&self) -> usize {
        cdim(&self.shell.l)
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

#[derive(Default)]
pub struct Basis {
    shells: Vec<BasisShell>,
}

impl Basis {
    pub fn new(atoms: &Vec<Atom>, shells: Vec<Vec<Shell>>) -> Self {
        let mut basis_shells: Vec<BasisShell> = Default::default();
        for atom in atoms {
            for shell in &shells[atom.z() as usize] {
                basis_shells.push(BasisShell::new(*atom.origin(), shell.clone()));
            }
        }

        Self {
            shells: basis_shells,
        }
    }

    /// return sperical harmonicas basis dimension
    pub fn sdim(self) -> usize {
        self.shells.iter().map(|s| s.sdim()).sum()
    }

    /// return sperical harmonicas basis dimension
    pub fn cdim(self) -> usize {
        self.shells.iter().map(|s| s.cdim()).sum()
    }

    // fn offset(&self, sn: &usize) -> usize {
    //     let mut offset = 0;
    //     for i in 0..*sn {
    //         offset += self.shells[i].dim();
    //     }
    //     offset
    // }
    //
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

        let l = i16::from(*lmax);

        let xyz = gaussian_layout(lmax);
        let mut mat = FMatrix::zero(usize::from(*lmax), usize::from(*lmax));

        for c in 0..cdim(lmax) {
            lx = xyz[c][0];
            ly = xyz[c][1];
            lz = xyz[c][2];

            for m in -l..=l {
                let j = i16::from(lx) + i16::from(ly) - m.abs();
                if j >= 0 && j % 2 == 0 {
                    // todo: body
                } else {
                    mat[((m + l) as usize, c)] = 0.0;
                }
            }
        }

        mat
    }
}

#[derive(Clone)]
struct CartesianBasisFunction {
    pub origin: [f64; 3],
    pub ml: [u8; 3],

    pub exps: Vec<f64>,
    pub coefs: Vec<f64>,
    pub norm: Vec<f64>,
}

impl CartesianBasisFunction {
    fn normalize(&mut self) {
        let l = self.ml[0];
        let m = self.ml[1];
        let n = self.ml[2];
        let l_total: f64 = (l + m + n).into();

        // have to cast to signed integer to prevent overflow when
        // either l, m, or n is 0
        let f_ijk: f64 = ((2 * i32::from(l) - 1).factorial2()
            * (2 * i32::from(m) - 1).factorial2()
            * (2 * i32::from(n) - 1).factorial2()) as f64;

        for e in &self.exps {
            let norm_e =
                2.0_f64.powf(2.0 * l_total + 1.5) * e.powf(l_total + 1.5) / f_ijk / PI.powf(1.5);
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

// Matrix cartesian_spherical_transformation(in int lmax) {
// int l, m, l2;
// // cartesian components
// int lx, ly, lz;
// int c, i, j, k, exponent;
// double s, s1, s2;
//
//
// l=lmax;
// auto xyz = GaussianLayout( lmax );
// auto matrix = slice!double( ldim( l ), cdim( l ) );
// {
// l2 = l*l;
// for(c=0; c<cdim(l); c++) {
// lx = xyz[c][0];
// ly = xyz[c][1];
// lz = xyz[c][2];
//
// for(m=-l; m<=l; m++) {
// j = lx+ly-abs(m);
// if (j>=0 && j%2==0) {
// j = j/2;
// s1 = 0.0;
//
// for(i=0; i<=(l-abs(m))/2; i++) {
// s2 = 0.0;
// for(k=0; k<=j; k++) {
// if( (m<0 && abs(abs(m)-lx)%2==1) ||
// (m>0 && abs(abs(m)-lx)%2==0) ) {
// exponent = (abs(m)-lx+2*k)/2;
// s = pow(-1.0, exponent)*sqrt(2.0);
// }
// else if (m==0 && lx%2==0) {
// exponent = k-lx/2;
// s = pow(-1.0, exponent);
// }
// else {
// s = 0.0;
// }
// s2 = s2 + n_over_k(j, k)*
// n_over_k(abs(m),(lx-2*k))*s;
// }
// s1 = s1 + n_over_k(l,i)*n_over_k(i,j)*
// pow(-1.0,i)*fact(2*l-2*i)/fact(l-abs(m)-2*i)*s2;
// }
// matrix[ m+l, c] =
// sqrt(1.0*(fact(2*lx)*fact(2*ly)*fact(2*lz)*fact(l)*fact(l-abs(m))) /
// (fact(lx)*fact(ly)*fact(lz)*fact(2*l)*fact(l+abs(m))))*
// s1 / (pow(2.0, l)*fact(l));
// }
// else
// matrix[ m+l, c] = 0.0;
// }
// }
// }
// return matrix;
// }

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
