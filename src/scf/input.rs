use libferric::HFType;

pub struct SCFInput {
    // HF type
    pub hf: HFType,

    // Convergence threshold
    pub e_threshold: f64,
    pub rms_threshold: f64,

    // iterations
    pub max_iter: usize,

    // diis
    pub diis_iter_start: usize,
    pub diis_dim_max: usize,
}

impl Default for SCFInput {
    fn default() -> Self {
        Self {
            hf: HFType::RHF,

            e_threshold: 1e-6,
            rms_threshold: 1e-12,

            max_iter: 40,

            diis_iter_start: 1,
            diis_dim_max: 6,
        }
    }
}
