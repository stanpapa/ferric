use libferric::HFType;

pub struct SCFInput {
    // HF type
    hf: HFType,

    // Convergence threshold
    conv_thresh: f64,
}
