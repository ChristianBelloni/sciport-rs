use ndarray::{concatenate, Array1, ArrayView1, Axis};
use num::{Complex, Float};

pub(crate) mod complex;

pub fn convolve<T: Float>(
    data: ArrayView1<Complex<T>>,
    window: ArrayView1<Complex<T>>,
) -> Array1<Complex<T>> {
    if window.len() > data.len() {
        return convolve(window, data);
    }
    let data = concatenate![Axis(0), Array1::zeros(window.len() - 1), data,];
    let mut w = window.view();
    w.invert_axis(Axis(0));

    data.windows(w.len())
        .into_iter()
        .map(|x| (&x * &w).sum())
        .collect()
}
