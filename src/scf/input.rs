use std::str::FromStr;

use libferric::HFType;
use serde_yaml::Value;

#[derive(Clone)]
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

            diis_iter_start: 2,
            diis_dim_max: 6,
        }
    }
}

impl SCFInput {
    pub fn parse(input: &Value) -> SCFInput {
        // initialise SCFInput with default values
        let mut scf = SCFInput::default();

        // overwrite defaults with input values
        for (key, value) in input.as_mapping().unwrap() {
            match key.as_str().unwrap().to_lowercase().as_str() {
                "hf" => scf.hf = HFType::from_str(value.as_str().unwrap()).unwrap(),
                "maxiter" => scf.max_iter = value.as_u64().unwrap() as usize,
                "diisiterstart" => scf.diis_iter_start = value.as_u64().unwrap() as usize,
                "diisdimmax" => scf.diis_dim_max = value.as_u64().unwrap() as usize,
                "thresholde" => scf.e_threshold = value.as_f64().unwrap(),
                "thresholdrms" => scf.rms_threshold = value.as_f64().unwrap(),
                _ => panic!("Unknown option: {:?}", key),
            }
        }

        // return SCFInput
        scf
    }
}
