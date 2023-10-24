use ndarray::{Array, Dimension};
use num::{Float, NumCast};

pub fn sinc<T: Float, D: Dimension>(x: Array<T, D>) -> Array<T, D> {
    let x = x.mapv(|a| <f64 as NumCast>::from(a).unwrap());
    x.mapv(|a| {
        if a == 0.0 {
            T::one()
        } else {
            T::from((a * std::f64::consts::PI).sin() / (a * std::f64::consts::PI)).unwrap()
        }
    })
}
