use ndarray::{array, Array1};
use numpy::Complex64;

/// Copy-pasted from scipy, this can probably be optimized
pub(crate) fn c_filt(
    mut b: Array1<Complex64>,
    mut a: Array1<Complex64>,
    signal: Array1<f64>,
    mut filter_state: Array1<Complex64>,
) -> Array1<f64> {
    let a0 = a[0];

    b = b / a0;
    a = a / a0;

    let mut ret = Array1::<f64>::zeros(signal.raw_dim());

    let len_signal = signal.len();

    for i in 0..len_signal {
        if b.len() > 1 {
            ret[i] = (filter_state[0] + b[0] * signal[i]).norm();
            for k in 1..(b.len() - 1) {
                filter_state[k] = filter_state[k + 1] + signal[i] * b[k] - ret[i] * a[k];
            }
            filter_state[b.len() - 1] = signal[i] * b[b.len() - 1] - ret[i] * a[a.len() - 1];
        } else {
            ret[i] = (signal[i] * b[0]).norm();
        }
    }
    ret
}

pub fn linear_filter(
    b: Array1<Complex64>,
    a: Array1<Complex64>,
    signal: Array1<f64>,
    filter_state: Array1<Complex64>,
) -> Array1<f64> {
    c_filt(b, a, signal, filter_state)
}
