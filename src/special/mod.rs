mod kv;
mod trig;
pub use kv::*;
use ndarray::Array1;
use num::Complex;
pub use trig::*;

pub fn i0(x: Array1<f64>) -> Result<Array1<Complex<f64>>, i32> {
    // TODO implement this using cephes
    x.mapv(|v| complex_bessel_rs::bessel_i::bessel_i(0.0, v.into()))
        .into_iter()
        .collect()
}
