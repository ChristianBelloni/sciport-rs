use ndarray::Array1;

use crate::signal::fir_filter_design::windows::{extend, len_guards, truncate};

pub fn barthann(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    let sym = sym.into().unwrap_or(true);
    fn _barthann(m: u64, sym: bool) -> Array1<f64> {
        if len_guards(m) {
            return Array1::zeros(m as usize);
        }

        let (m, needs_trunc) = extend(m, sym);

        let n = (0..m)
            .into_iter()
            .map(|a| a as f64)
            .collect::<Array1<f64>>();
        let m = m as f64;
        let fac = (n / (m - 1.0) - 0.5).mapv(|a| a.abs());
        let w =
            0.62 - 0.48 * fac.clone() + 0.38 * (2.0 * std::f64::consts::PI * fac).mapv(|a| a.cos());
        truncate(w, needs_trunc)
    }
    _barthann(m, sym)
}
