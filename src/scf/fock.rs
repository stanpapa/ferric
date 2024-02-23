use libferric::linear_algebra::{
    matrix::FMatrix, matrix_container::FMatrixSymContainer, traits::Dot,
};

/// build Fock matrix as
///     RHF: Fμν = Hμν + \sum_{ρσ} Dρσ [ (μν|ρσ) - 0.5 (μρ|νσ) ]
///     UHF: Fμν = Hμν + \sum_{ρσ} Dρσ [ 2 (μν|ρσ) - (μρ|νσ) ]
pub fn fock(d: &FMatrix, h: &FMatrix, eri: &FMatrixSymContainer, facj: f64, facx: f64) -> FMatrix {
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

    h.clone() + facj * j - facx * k
}
