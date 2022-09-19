use ferric_lib::HFType;

pub struct SCFInput {
    // HF type
    hf: HFType,

    // Convergence threshold
    conv_thresh: f64,
}
