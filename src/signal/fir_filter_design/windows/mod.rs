mod barthann;
mod utils;

use ndarray::{array, s, Array1};
pub(crate) use utils::*;

pub use barthann::*;

pub enum WindowType {
    Boxcar,
    Triang,
    Blackman,
    Hamming,
    Hann,
    Bartlett,
    FlatTop,
    Parzen,
    Bohman,
    BlackmanHarris,
    Nuttall,
    Barthann,
    Cosine,
    Exponential,
    Tukey,
    Taylor,
    Lanczos,
    Kaiser { beta: f64 },
    KaiserBesselDerived { beta: f64 },
    Gaussian { std_dev: f64 },
    GeneralCosine { coeffs: Array1<f64> },
    GeneralGaussian { power: f64, width: f64 },
    Dpss { half_bandwidth: f64 },
    Chebwin { attenuation: f64 },
}

pub fn get_window(window: WindowType, nx: u64, fftbins: impl Into<Option<bool>>) {
    let fftbins = fftbins.into().unwrap_or(true);
    let m = nx;
    let sym = !fftbins;
    match window {
        WindowType::Boxcar => boxcar(m, sym),
        WindowType::Triang => triang(m, sym),
        WindowType::Blackman => blackman(m, sym),
        WindowType::Hamming => hamming(m, sym),
        WindowType::Hann => hann(m, sym),
        WindowType::Bartlett => bartlett(m, sym),
        WindowType::FlatTop => flattop(m, sym),
        WindowType::Parzen => parzen(m, sym),
        WindowType::Bohman => todo!(),
        WindowType::BlackmanHarris => todo!(),
        WindowType::Nuttall => todo!(),
        WindowType::Barthann => barthann(nx, fftbins),
        WindowType::Cosine => todo!(),
        WindowType::Exponential => todo!(),
        WindowType::Tukey => todo!(),
        WindowType::Taylor => todo!(),
        WindowType::Lanczos => todo!(),
        WindowType::Kaiser { beta } => todo!(),
        WindowType::KaiserBesselDerived { beta } => todo!(),
        WindowType::Gaussian { std_dev } => todo!(),
        WindowType::GeneralCosine { coeffs } => todo!(),
        WindowType::GeneralGaussian { power, width } => todo!(),
        WindowType::Dpss { half_bandwidth } => todo!(),
        WindowType::Chebwin { attenuation } => todo!(),
    };
}

pub fn boxcar(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _boxcar(m: u64, sym: bool) -> Array1<f64> {
        if len_guards(m) {
            return Array1::zeros(m as usize);
        }

        let (m, needs_trunc) = extend(m, sym);

        let w = Array1::ones(m as usize);
        truncate(w, needs_trunc)
    }

    _boxcar(m, sym.into().unwrap_or(true))
}

pub fn triang(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _triang(m: u64, sym: bool) -> Array1<f64> {
        if len_guards(m) {
            return Array1::ones(m as usize);
        }

        let (m, needs_trunc) = extend(m, sym);

        let n = 1..((m + 1) / 2 + 1);
        let n = n.map(|a| a as _).collect::<Array1<f64>>();

        let w = if m % 2 == 0 {
            let w = (2.0 * n - 1.0) / (m as f64);
            let reversed = w.slice(s![..;-1]);
            let mut w = w.to_vec();
            let reversed = reversed.to_vec();
            w.extend(&reversed);
            Array1::from(w)
        } else {
            let w = 2.0 * n / ((m as f64) + 1.0);
            let reversed = w.slice(s![-2..;-2]);
            let mut w = w.to_vec();
            let reversed = reversed.to_vec();
            w.extend(&reversed);
            Array1::from(w)
        };

        truncate(w, needs_trunc)
    }

    _triang(m, sym.into().unwrap_or(true))
}

pub fn general_cosine(m: u64, a: Array1<f64>, sym: impl Into<Option<bool>>) -> Array1<f64> {
    let sym = sym.into().unwrap_or(true);
    fn _general_cosine(m: u64, a: Array1<f64>, sym: bool) -> Array1<f64> {
        if len_guards(m) {
            return Array1::zeros(m as usize);
        }

        let (m, needs_trunc) = extend(m, sym);

        let fac = Array1::linspace(-std::f64::consts::PI, std::f64::consts::PI, m as _);
        let mut w = Array1::<f64>::zeros((m as _,));

        for (k, a) in a.iter().enumerate() {
            w = w + fac.map(|e| (e * (k as f64)).cos()) * (*a);
        }
        truncate(w, needs_trunc)
    }
    _general_cosine(m, a, sym)
}

pub fn blackman(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    let sym = sym.into().unwrap_or(true);

    fn _blackman(m: u64, sym: bool) -> Array1<f64> {
        general_cosine(m, array![0.42, 0.50, 0.08], sym)
    }

    _blackman(m, sym)
}

pub fn general_hamming(m: u64, alpha: f64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _general_hamming(m: u64, alpha: f64, sym: bool) -> Array1<f64> {
        general_cosine(m, array![alpha, 1.0 - alpha], sym)
    }

    _general_hamming(m, alpha, sym.into().unwrap_or(true))
}

pub fn hamming(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _hamming(m: u64, sym: bool) -> Array1<f64> {
        general_hamming(m, 0.54, sym)
    }
    _hamming(m, sym.into().unwrap_or(true))
}

pub fn hann(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _hann(m: u64, sym: bool) -> Array1<f64> {
        general_hamming(m, 0.5, sym)
    }
    _hann(m, sym.into().unwrap_or(true))
}

pub fn bartlett(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _bartlett(m: u64, sym: bool) -> Array1<f64> {
        if len_guards(m) {
            return Array1::ones((m as usize,));
        }

        let (m, needs_trunc) = extend(m, sym);

        let n = Array1::range(0.0, m as f64, 1.0);

        let filter_map = n.map(|&a| a <= (m as f64 - 1.0) / 2.0);
        let left = 2.0 * n.clone() / (m as f64 - 1.0);
        let right = 2.0 - 2.0 * n / (m as f64 - 1.0);
        let mut w = Array1::<f64>::zeros((m as usize,));
        for (i, v) in filter_map.into_iter().enumerate() {
            w[i] = if v { left[i] } else { right[i] };
        }
        truncate(w, needs_trunc)
    }
    _bartlett(m, sym.into().unwrap_or(true))
}

pub fn flattop(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _flattop(m: u64, sym: bool) -> Array1<f64> {
        let a = array![
            0.21557895,
            0.41663158,
            0.277263158,
            0.083578947,
            0.006947368
        ];
        general_cosine(m, a, sym)
    }

    _flattop(m, sym.into().unwrap_or(true))
}

fn extract<T: Clone>(condition: Array1<bool>, arr: Array1<T>) -> Array1<T> {
    let size = condition.iter().filter(|&a| *a).count();
    let mut result = Vec::with_capacity(size);
    let mut arr = arr.to_vec();
    let mut removed = 0;
    for (i, v) in condition.into_iter().enumerate() {
        if v {
            result.push(arr.remove(i - removed));
            removed += 1;
        }
    }
    result.into()
}

pub fn parzen(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _parzen(m: u64, sym: bool) -> Array1<f64> {
        if len_guards(m) {
            return Array1::ones((m as usize,));
        }

        let (m, needs_trunc) = extend(m, sym);
        let m = m as f64;
        let n = Array1::range(-(m - 1.0) / 2.0, (m - 1.0) / 2.0 + 0.5, 1.0);

        let na = extract(n.map(|&v| v < -(m - 1.0) / 4.0), n.clone());
        let nb = extract(n.map(|&v| v.abs() <= (m - 1.0) / 4.0), n);

        let wa = (2.0 * (1.0 - na.mapv(f64::abs) / (m / 2.0))).mapv(|v| v.powi(3));

        let wb = (1.0 - 6.0 * (nb.mapv(f64::abs) / (m / 2.0))).mapv(|v| v.powi(2))
            + (6.0 * (nb.mapv(f64::abs) / (m / 2.0))).mapv(|v| v.powi(3));

        let mut w = wa.clone().to_vec();
        w.extend(&wb);
        w.extend(wa.slice(s![..;-1]));

        truncate(w.into(), needs_trunc)
    }
    _parzen(m, sym.into().unwrap_or(true))
}
