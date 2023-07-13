use ndarray::{Array, Array0, Axis, Dimension};
use num::Num;

pub fn convolve<N: Num + Copy, D>(
    in1: Array<N, D>,
    in2: Array<N, D>,
    method: ConvolveMethod,
    mode: ConvolveMode,
) {
}

pub fn fft_convolve<N: Num + Copy + Default, D: Dimension>(
    in1: Array<N, D>,
    in2: Array<N, D>,
    mode: ConvolveMode,
) -> Array<N, D> {
    if in1.ndim() == 0 {
        return in1 * in2;
    }

    if in1.len() == 0 || in2.len() == 0 {
        return Default::default();
    }

    todo!()
}

fn _init_freq_conv_axes<N: Num + Copy + Default, D: Dimension>(
    in1: Array<N, D>,
    in2: Array<N, D>,
    mode: ConvolveMode,
    axes: impl Into<Option<Vec<Axis>>>,
    sorted_axes: bool,
) {
    let axes: Option<_> = axes.into();
    let s1 = in1.shape();
    let s2 = in2.shape();

    let noaxes = axes.is_none();
}

pub enum ConvolveMethod {
    Auto,
    Fft,
    Direct,
}

pub enum ConvolveMode {
    Full,
    Valid,
    Same,
}
