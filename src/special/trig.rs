use ndarray::{Array, Dimension};

pub fn sinc<D: Dimension>(x: Array<f64, D>) -> Array<f64, D> {
    x.mapv(|a| {
        if a == 0.0 {
            1.0
        } else {
            (a * std::f64::consts::PI).sin() / (a * std::f64::consts::PI)
        }
    })
}
