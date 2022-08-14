// use super::gto_bases::*;
use crate::constants::PI;
use crate::math::functions::factorial2;

#[derive(Clone)]
pub struct Shell {
    l: i32,
    exps: Vec<f64>,
    coefs: Vec<f64>,
}

impl Default for Shell {
    fn default() -> Self {
        Shell {
            l: -1,
            exps: vec![],
            coefs: vec![],
        }
    }
}

impl Shell {
    pub fn new(l: i32, exps: Vec<f64>, coefs: Vec<f64>) -> Shell {
        Shell {
            l: l,
            exps: exps,
            coefs: coefs,
        }
    }
}

pub struct BasisShell {
    origin: [f64; 3],
    shell: Shell,

    cbf: Vec<CartesianBasisFunction>,
}

impl BasisShell {
    pub fn new(origin: [f64; 3], shell: Shell) -> BasisShell {
        let mut basis_shell = BasisShell {
            origin: origin,
            shell: shell,
            cbf: vec![],
        };
        basis_shell.generate_cartesian_shell();
        return basis_shell;
    } // new

    pub fn dim(&self) -> i32 {
        return 2 * self.shell.l + 1;
    }
    pub fn cdim(&self) -> i32 {
        return (self.shell.l + 1) * (self.shell.l + 2) / 2;
    }

    fn generate_cartesian_shell(&mut self) {
        let gaussians: Vec<[i32; 3]> = gaussian_components(&self.shell.l);
        self.cbf = Vec::with_capacity(gaussians.len());
        for g in gaussians.iter() {
            self.cbf.push(CartesianBasisFunction {
                origin: self.origin,
                ml: *g,
                exps: (*self.shell.exps).to_vec(),
                coefs: (*self.shell.coefs).to_vec(),
                norm: vec![],
            });
        } // g
        for i in 0..self.cbf.len() {
            self.cbf[i].normalize();
        } // f
    } // fn generate_cartesian_shell
} // impl BasisShell

pub struct Basis {
    shells: Vec<BasisShell>,
}

impl Basis {
    pub fn dim(self) -> i32 {
        let mut dim = 0;
        for shell in self.shells {
            dim += shell.dim();
        }
        return dim;
    } // fn dim

    pub fn cdim(self) -> i32 {
        let mut dim = 0;
        for shell in self.shells {
            dim += shell.cdim();
        }
        return dim;
    } // fn cdim
} // impl Basis

struct CartesianBasisFunction {
    pub origin: [f64; 3],
    pub ml: [i32; 3],

    pub exps: Vec<f64>,
    pub coefs: Vec<f64>,
    pub norm: Vec<f64>,
}

impl CartesianBasisFunction {
    fn normalize(&mut self) {
        // fn normalize(mut self) {
        let i = self.ml[0];
        let j = self.ml[1];
        let k = self.ml[2];
        let l: f64 = (i + j + k).into();
        self.norm = Vec::with_capacity(self.exps.len());
        for a in self.exps.iter() {
            let Fijk: f64 =
                (factorial2(2 * i - 1) * factorial2(2 * j - 1) * factorial2(2 * k - 1)).into();
            self.norm.push(
                (2.0_f64.powf(2.0 * l + 1.5) * a.powf(l + 1.5) / (PI.powf(1.5) * Fijk)).sqrt(),
            );
        } // a

        // TODO prefactor?

        let mut N: f64 = 0.0;
        for bra in 0..self.exps.len() {
            for ket in 0..self.exps.len() {
                N += self.norm[bra] * self.norm[ket] * self.coefs[bra] * self.coefs[ket] / 1.0;
            } // ket
        } // bra

        N = N.powf(-0.5);
        for c in 0..self.coefs.len() {
            self.coefs[c] *= N;
        } // c
    } // normalize
} // impl CartesianBasisFunction

fn gaussian_components(l: &i32) -> Vec<[i32; 3]> {
    let n = (l + 1) * (l + 2) / 2;
    // l = &2;
    let mut components: Vec<[i32; 3]> = Vec::with_capacity(n as usize);
    let mut i = *l;
    let mut j = 0;
    let mut k = 0;
    let mut count = 0;
    while {
        components.push([i, j, k]);
        if k < l - i {
            // j = j - 1;
            j -= 1;
            k += 1;
        } else {
            i -= 1;
            j = l - i;
            k = 0;
        }
        count += 1;
        count < n
    } {}

    return components;
} // fn gaussian_components

// #[cfg(test)]
// mod tests {
//   use super::*;

//   #[test]
//   fn test_CBF_normalize() {
//     let mut cbf = CartesianBasisFunction {
//       origin: [0.0, 0.0, 0.0],
//       ml: [1,0,0],
//       exps: vec![4.2, 2.1, 1.0],
//       coefs: vec![1.0, 1.0, 1.0],
//       norm: vec![0.0, 0.0, 0.0],
//     };
//     cbf.normalize();
//   }
// }
