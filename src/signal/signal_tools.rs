use ndarray::Array1;
use num::{Complex, Float};

/// Copy-pasted from scipy, this can probably be optimized
pub(crate) fn c_filt<T: Float>(
    mut b: Array1<Complex<T>>,
    mut a: Array1<Complex<T>>,
    signal: Array1<Complex<T>>,
    mut filter_state: Array1<Complex<T>>,
) -> Array1<Complex<T>> {
    let a0 = a[0];

    b = b.mapv(|a| a / a0);
    a = a.mapv(|a| a / a0);

    let mut ret = Array1::<Complex<T>>::zeros(signal.raw_dim());

    let len_signal = signal.len();

    for i in 0..len_signal {
        if b.len() > 1 {
            ret[i] = filter_state[0] + signal[i] * b[0];
            for n in 1..(b.len() - 1) {
                filter_state[n - 1] = filter_state[n] + signal[i] * b[n] - ret[i] * a[n];
            }
            filter_state[b.len() - 2] = signal[i] * b[b.len() - 1] - ret[i] * a[b.len() - 1];
        }
    }
    ret
}

pub fn linear_filter<T: Float>(
    b: Array1<Complex<T>>,
    a: Array1<Complex<T>>,
    signal: Array1<Complex<T>>,
    filter_state: Array1<Complex<T>>,
) -> Array1<Complex<T>> {
    c_filt(b, a, signal, filter_state)
}
