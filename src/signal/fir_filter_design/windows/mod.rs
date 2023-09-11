mod barthann;
mod utils;

pub use barthann::*;
use ndarray::{array, s, Array1};
use std::f64::consts::PI;
pub(crate) use utils::*;

use crate::if_len_guard;

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
    Exponential {
        center: Option<f64>,
        tau: Option<f64>,
    },
    Tukey {
        alpha: Option<f64>,
    },
    Taylor,
    Lanczos,
    Kaiser {
        beta: f64,
    },
    KaiserBesselDerived {
        beta: f64,
    },
    Gaussian {
        std_dev: f64,
    },
    GeneralCosine {
        coeffs: Array1<f64>,
    },
    GeneralGaussian {
        power: f64,
        width: f64,
    },
    Dpss {
        half_bandwidth: f64,
    },
    Chebwin {
        attenuation: f64,
    },
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
        WindowType::Bohman => bohman(m, sym),
        WindowType::BlackmanHarris => blackmanharris(m, sym),
        WindowType::Nuttall => nuttall(m, sym),
        WindowType::Barthann => barthann(nx, fftbins),
        WindowType::Cosine => cosine(m, sym),
        WindowType::Exponential { center, tau } => exponential(m, center, tau, sym),
        WindowType::Tukey { alpha } => tukey(m, alpha, sym),
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
        if_len_guard!(m);

        let (m, needs_trunc) = extend(m, sym);

        let w = Array1::ones(m as usize);
        truncate(w, needs_trunc)
    }

    _boxcar(m, sym.into().unwrap_or(true))
}

pub fn triang(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _triang(m: u64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);

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

pub fn general_cosine(
    m: u64,
    a: impl Into<Array1<f64>>,
    sym: impl Into<Option<bool>>,
) -> Array1<f64> {
    let sym = sym.into().unwrap_or(true);
    fn _general_cosine(m: u64, a: Array1<f64>, sym: bool) -> Array1<f64> {
        if_len_guard!(m);

        let (m, needs_trunc) = extend(m, sym);

        let fac = Array1::linspace(-PI, PI, m as _);
        let mut w = Array1::<f64>::zeros((m as _,));

        for (k, a) in a.iter().enumerate() {
            w = w + fac.map(|e| (e * (k as f64)).cos()) * (*a);
        }
        truncate(w, needs_trunc)
    }
    _general_cosine(m, a.into(), sym)
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
        if_len_guard!(m);

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
        if_len_guard!(m);

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

        truncate(w, needs_trunc)
    }
    _parzen(m, sym.into().unwrap_or(true))
}

pub fn bohman(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _bohman(m: u64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);

        let (m, needs_trunc) = extend(m, sym);

        let fac = Array1::linspace(-1.0, 1.0, m as _)
            .slice(s![1..-1])
            .mapv(f64::abs);
        let temp = (1.0 - fac.clone()) * fac.mapv(|a| (PI * a).cos())
            + 1.0 / PI * fac.mapv(|a| (PI * a).sin());
        let mut w = Vec::with_capacity(temp.len());
        w.push(0.0);
        w.extend(temp);
        w.push(0.0);
        truncate(w, needs_trunc)
    }
    _bohman(m, sym.into().unwrap_or(true))
}

pub fn blackmanharris(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _blackmanharris(m: u64, sym: bool) -> Array1<f64> {
        general_cosine(m, vec![0.35875, 0.48829, 0.14128, 0.01168], sym)
    }

    _blackmanharris(m, sym.into().unwrap_or(true))
}

pub fn nuttall(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _nuttall(m: u64, sym: bool) -> Array1<f64> {
        general_cosine(m, vec![0.3635819, 0.4891775, 0.1365995, 0.0106411], sym)
    }

    _nuttall(m, sym.into().unwrap_or(true))
}

pub fn cosine(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _cosine(m: u64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);
        let (m, needs_trunc) = extend(m, sym);

        let w = Array1::range(0.0, m as _, 1.0);
        let m = m as f64;
        let w = (PI / m * w + 0.5).mapv(f64::sin);

        truncate(w, needs_trunc)
    }
    _cosine(m, sym.into().unwrap_or(true))
}

pub fn exponential(
    m: u64,
    center: impl Into<Option<f64>>,
    tau: impl Into<Option<f64>>,
    sym: impl Into<Option<bool>>,
) -> Array1<f64> {
    fn _exponential(m: u64, center: Option<f64>, tau: f64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);
        let (m, needs_trunc) = extend(m, sym);

        let center = if let Some(center) = center {
            center
        } else {
            (m as f64) - 1.0 / 2.0
        };

        let n = Array1::range(0.0, m as f64, 1.0);
        let w = (-(n - center).mapv(f64::abs) / tau).mapv(f64::exp);

        truncate(w, needs_trunc)
    }

    _exponential(
        m,
        center.into(),
        tau.into().unwrap_or(1.0),
        sym.into().unwrap_or(true),
    )
}

pub fn tukey(m: u64, alpha: impl Into<Option<f64>>, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _tukey(m: u64, alpha: f64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);

        if alpha <= 0.0 {
            return Array1::ones((m as usize,));
        } else if alpha >= 1.0 {
            return hann(m, sym);
        }

        let (m, needs_trunc) = extend(m, sym);
        let m = m as f64;
        let n = Array1::range(0.0, m, 1.0);

        let width = (alpha * (m - 1.0) / 2.0).floor() as usize;

        let n1 = n.slice(s![0..(width + 1)]);
        let n2 = n.slice(s![(width + 1)..((m as usize) - width - 1)]);
        let n3 = n.slice(s![(m as usize - width - 1)..]);

        let w1 =
            0.5 * (1.0 + (PI * (-1.0 + 2.0 * n1.to_owned() / alpha / (m - 1.0))).mapv(f64::cos));
        let w2 = Array1::<f64>::zeros((n2.len(),));
        let w3 = 0.5
            * (1.0
                + (PI * (-2.0 / alpha + 1.0 + 2.0 * n3.to_owned() / alpha / (m - 1.0)))
                    .mapv(f64::cos));

        let mut w = Vec::with_capacity(w1.len() + w2.len() + w3.len());
        w.extend(w1);
        w.extend(w2);
        w.extend(w3);
        truncate(w, needs_trunc)
    }

    _tukey(m, alpha.into().unwrap_or(0.5), sym.into().unwrap_or(true))
}
