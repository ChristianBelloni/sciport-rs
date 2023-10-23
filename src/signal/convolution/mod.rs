use ndarray::{Array, Axis, Dimension};
use num::Num;

#[allow(unused)]
pub fn convolve<N: Num + Copy, D>(
    in1: Array<N, D>,
    in2: Array<N, D>,
    method: ConvolveMethod,
    mode: ConvolveMode,
) {
}

#[allow(unused)]
pub fn fft_convolve<N: Num + Copy + Default, D: Dimension>(
    in1: Array<N, D>,
    in2: Array<N, D>,
    mode: ConvolveMode,
) -> Array<N, D> {
    if in1.ndim() == 0 {
        return in1 * in2;
    }

    if in1.is_empty() || in2.is_empty() {
        return Default::default();
    }

    todo!()
}

#[allow(unused)]
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

#[allow(unused)]
pub enum ConvolveMethod {
    Auto,
    Fft,
    Direct,
}

#[allow(unused)]
pub enum ConvolveMode {
    Full,
    Valid,
    Same,
}
