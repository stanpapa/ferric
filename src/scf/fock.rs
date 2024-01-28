use ferric_lib::linear_algebra::{
    matrix::FMatrix, matrix_container::FMatrixSymContainer, traits::Dot,
};

// TODO: FIX UHF FOCK BUILD

/// build Fock matrix as
///     RHF: Fμν = Hμν + \sum_{ρσ} Dρσ [ (μν|ρσ) - 0.5 (μρ|νσ) ]
///     UHF: Fμν = Hμν + \sum_{ρσ} Dρσ [ (μν|ρσ) - (μρ|νσ) ]
pub fn fock(h: &FMatrix, d: &[FMatrix; 2], eri: &FMatrixSymContainer) -> [FMatrix; 2] {
    let uhf = d[1].rows > 0;

    let mut f = [FMatrix::default(), FMatrix::default()];
    if uhf {
        f = fock_uhf(h, d, eri);
    } else {
        f[0] = fock_rhf(h, &d[0], eri);
    }

    f
}

/// build RHF Fock matrix as
///     Fμν = Hμν + \sum_{ρσ} Dρσ [ (μν|ρσ) - 0.5 (μρ|νσ) ]
fn fock_rhf(h: &FMatrix, d: &FMatrix, eri: &FMatrixSymContainer) -> FMatrix {
    let dim = h.rows;

    let mut j = FMatrix::zero(dim, dim);
    for μ in 0..dim {
        for ν in 0..dim {
            // only lower-triangle is stored
            let μ1 = std::cmp::max(μ, ν);
            let ν1 = std::cmp::min(μ, ν);
            j[(μ, ν)] = d.dot(&eri[(μ1, ν1)]);
        }
    }

    let mut k = FMatrix::zero(dim, dim);
    for μ in 0..dim {
        for ν in 0..dim {
            for ρ in 0..dim {
                // only lower-triangle is stored
                let μ1 = std::cmp::max(μ, ρ);
                let ρ1 = std::cmp::min(μ, ρ);
                for σ in 0..dim {
                    k[(μ, ν)] += d[(ρ, σ)] * eri[(μ1, ρ1)][(ν, σ)];
                }
            }
        }
    }

    h.clone() + j - 0.5 * k
}

/// build Fock matrix as
///     UHF: Fμν = Hμν + \sum_{ρσ} Dρσ [ (μν|ρσ) - (μρ|νσ) ]
fn fock_uhf(h: &FMatrix, d: &[FMatrix; 2], eri: &FMatrixSymContainer) -> [FMatrix; 2] {
    let dim = h.rows;

    // -----------------------
    // alpha
    // -----------------------
    let mut j_a = FMatrix::zero(dim, dim);
    for μ in 0..dim {
        for ν in 0..dim {
            // only lower-triangle is stored
            let μ1 = std::cmp::max(μ, ν);
            let ν1 = std::cmp::min(μ, ν);
            j_a[(μ, ν)] = d[0].dot(&eri[(μ1, ν1)]);
        }
    }

    let mut k_a = FMatrix::zero(dim, dim);
    for μ in 0..dim {
        for ν in 0..dim {
            for ρ in 0..dim {
                // only lower-triangle is stored
                let μ1 = std::cmp::max(μ, ρ);
                let ρ1 = std::cmp::min(μ, ρ);
                for σ in 0..dim {
                    k_a[(μ, ν)] += d[0][(ρ, σ)] * eri[(μ1, ρ1)][(ν, σ)];
                }
            }
        }
    }

    let f_a = h.clone() + j_a - k_a;

    // -----------------------
    // beta
    // -----------------------
    let mut j_b = FMatrix::zero(dim, dim);
    for μ in 0..dim {
        for ν in 0..dim {
            // only lower-triangle is stored
            let μ1 = std::cmp::max(μ, ν);
            let ν1 = std::cmp::min(μ, ν);
            j_b[(μ, ν)] = d[1].dot(&eri[(μ1, ν1)]);
        }
    }

    let mut k_b = FMatrix::zero(dim, dim);
    for μ in 0..dim {
        for ν in 0..dim {
            for ρ in 0..dim {
                // only lower-triangle is stored
                let μ1 = std::cmp::max(μ, ρ);
                let ρ1 = std::cmp::min(μ, ρ);
                for σ in 0..dim {
                    k_b[(μ, ν)] += d[1][(ρ, σ)] * eri[(μ1, ρ1)][(ν, σ)];
                }
            }
        }
    }

    let f_b = h.clone() + j_b - k_b;

    [f_a, f_b]
}
