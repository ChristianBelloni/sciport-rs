use ndarray::{Array, Axis, Dimension};
use num::Num;

#[allow(unused)]
pub fn convolve<N: Num + Copy + Default, D: Dimension>(
    mut in1: Array<N, D>,
    mut in2: Array<N, D>,
    method: ConvolveMethod,
    mode: ConvolveMode,
) -> Array<N, D> {
    if inputs_swap_needed(mode, in1.shape(), in2.shape()) {
        std::mem::swap(&mut in1, &mut in2);
    }

    match method {
        ConvolveMethod::Auto => todo!(),
        ConvolveMethod::Fft => fft_convolve(in1, in2, mode),
        ConvolveMethod::Direct => todo!(),
    }
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

fn inputs_swap_needed(mode: ConvolveMode, shape1: &[usize], shape2: &[usize]) -> bool {
    debug_assert_eq!(shape1.len(), shape2.len());

    if !matches!(mode, ConvolveMode::Valid) {
        return false;
    }

    shape1.iter().zip(shape2).all(|(a, b)| a >= b) | shape2.iter().zip(shape1).all(|(a, b)| a >= b)
}

#[allow(unused)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ConvolveMethod {
    Auto,
    Fft,
    Direct,
}

#[allow(unused)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ConvolveMode {
    Full,
    Valid,
    Same,
}
