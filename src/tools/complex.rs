use ndarray::{Array, Dimension};
use num::{Complex, Float};

pub fn normalize_zeros<D: Dimension, T: Float>(a: Array<Complex<T>, D>) -> Array<Complex<T>, D> {
    a.mapv(|mut a| {
        if a.im == T::neg_zero() {
            a.im = T::zero()
        }
        a
    })
}
