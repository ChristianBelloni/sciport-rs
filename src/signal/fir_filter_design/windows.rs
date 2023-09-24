use ndarray::Array1;
use ndarray::{array, s, NewAxis};
use std::{f64::consts::PI, ops::Div};

use crate::{if_len_guard, special::sinc};

#[allow(unused)]
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
    Taylor {
        nbar: Option<u64>,
        sll: Option<f64>,
        norm: Option<bool>,
    },
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

#[allow(unused)]
pub fn get_window(window: WindowType, nx: u64, fftbins: impl Into<Option<bool>>) -> Array1<f64> {
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
        WindowType::Taylor { nbar, sll, norm } => taylor(m, nbar, sll, norm, sym),
        WindowType::Lanczos => lanczos(m, sym),
        WindowType::Kaiser { beta } => todo!(),
        WindowType::KaiserBesselDerived { beta } => todo!(),
        WindowType::Gaussian { std_dev } => todo!(),
        WindowType::GeneralCosine { coeffs } => todo!(),
        WindowType::GeneralGaussian { power, width } => todo!(),
        WindowType::Dpss { half_bandwidth } => todo!(),
        WindowType::Chebwin { attenuation } => todo!(),
    }
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

        let n = Array1::range(1.0, (m as f64 + 1.0).div(2.0).floor() + 1.0, 1.0);

        let w = if m % 2 == 0 {
            let w = (2.0 * n - 1.0) / (m as f64);
            let mut reversed = w.to_vec();
            reversed.reverse();
            let mut w = w.to_vec();
            w.extend(&reversed);
            w
        } else {
            let w = 2.0 * n / ((m as f64) + 1.0);
            let mut reversed = w.to_vec();
            reversed.remove(reversed.len() - 1);
            reversed.reverse();
            let mut w = w.to_vec();
            w.extend(&reversed);
            w
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
        let wa = 2.0 * (1.0 - na.mapv(f64::abs) / (m / 2.0)).mapv(|v| v.powi(3));
        let wb = 1.0 - 6.0 * (nb.mapv(f64::abs) / (m / 2.0)).mapv(|v| v.powi(2))
            + 6.0 * (nb.mapv(f64::abs) / (m / 2.0)).mapv(|v| v.powi(3));
        let mut w = wa.clone().to_vec();
        w.extend(&wb);
        w.extend(wa.to_vec().into_iter().rev());

        truncate(w, needs_trunc)
    }
    _parzen(m, sym.into().unwrap_or(true))
}
#[allow(clippy::reversed_empty_ranges)]
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

pub fn barthann(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    let sym = sym.into().unwrap_or(true);
    fn _barthann(m: u64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);

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

pub fn cosine(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _cosine(m: u64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);
        let (m, needs_trunc) = extend(m, sym);

        let w = Array1::range(0.0, m as _, 1.0) + 0.5;
        let m = m as f64;
        let w = (PI / m * w).mapv(f64::sin);

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
            ((m as f64) - 1.0) / 2.0
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
        let w2 = Array1::<f64>::ones((n2.len(),));
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

pub fn taylor(
    m: u64,
    nbar: impl Into<Option<u64>>,
    sll: impl Into<Option<f64>>,
    norm: impl Into<Option<bool>>,
    sym: impl Into<Option<bool>>,
) -> Array1<f64> {
    fn _taylor(m: u64, nbar: u64, sll: f64, _norm: bool, sym: bool) -> Array1<f64> {
        if_len_guard!(m);
        let (m, needs_trunc) = extend(m, sym);

        let b = 10.0_f64.powf(sll / 20.0);
        let a = b.acosh() / PI;
        let nbar = nbar as f64;
        let s2 = nbar.powi(2) / (a.powi(2) + (nbar - 0.5).powi(2));
        let ma = Array1::range(1.0, nbar, 1.0);
        let mut fm = Array1::<f64>::zeros((nbar as usize - 1,));
        let signs =
            Array1::<f64>::from_shape_fn(ma.raw_dim(), |i| if i % 2 == 0 { 1.0 } else { -1.0 });
        let m2 = ma.clone() * ma.clone();

        for (mi, _) in ma.iter().enumerate() {
            let numer = signs[mi]
                * (1.0 - m2[mi] / s2 / (a.powi(2) + (ma.clone() - 0.5).mapv(|a| a.powi(2))))
                    .product();
            let denom = 2.0
                * (1.0 - m2[mi] / m2.slice(s![..mi]).to_owned()).product()
                * (1.0 - m2[mi] / m2.slice(s![mi + 1..]).to_owned()).product();
            fm[mi] = numer / denom;
        }

        let w_f = |n: Array1<f64>| {
            let rhs = 2.0 * PI * ma.slice(s![.., NewAxis]).to_owned() * (n - m as f64 / 2.0 + 0.5)
                / m as f64;
            let rhs = rhs.mapv(f64::cos);
            let ret = 1.0 + 2.0 * fm.dot(&rhs);
            ret
        };

        let w = w_f(Array1::range(0.0, m as f64, 1.0));

        truncate(w, needs_trunc)
    }

    _taylor(
        m,
        nbar.into().unwrap_or(4),
        sll.into().unwrap_or(30.0),
        norm.into().unwrap_or(true),
        sym.into().unwrap_or(true),
    )
}

pub fn lanczos(m: u64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _lanczos(m: u64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);
        let (m, needs_trunc) = extend(m, sym);

        fn _calc_right_side(n: u64, m: u64) -> Array1<f64> {
            let m = m as f64;
            let n = n as f64;
            sinc(2.0 * Array1::range(n, m, 1.0) / (m - 1.0) - 1.0)
        }
        let w = if m % 2 == 0 {
            let wh = _calc_right_side(m / 2, m);
            let mut lhs = wh.clone().to_vec();
            lhs.reverse();
            lhs.extend(wh);
            lhs
        } else {
            let wh = _calc_right_side((m + 1) / 2, m);
            let mut lhs = wh.clone().to_vec();
            lhs.reverse();
            lhs.push(1.0);
            lhs.extend(wh);
            lhs
        };

        truncate(w, needs_trunc)
    }
    _lanczos(m, sym.into().unwrap_or(true))
}

pub fn kaiser(m: u64, beta: f64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _kaiser(m: u64, beta: f64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);

        let (m, needs_trunc) = extend(m, sym);

        let n = Array1::range(0.0, m as f64, 1.0);
        let alpha = (m as f64 - 1.0) / 2.0;
        let l = beta * (1.0 - ((n - alpha) / alpha).mapv(|a| a.powi(2))).mapv(|a| a.sqrt());
        let l = crate::special::i0(l).unwrap();
        let r = crate::special::i0(array![beta]).unwrap();

        let w = l / r;
        let w = w.mapv(|a| a.norm());

        truncate(w, needs_trunc)
    }
    _kaiser(m, beta, sym.into().unwrap_or(true))
}

pub fn kaiser_bessel_derived(m: u64, beta: f64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _kaiser_bessel_derived(m: u64, beta: f64, _sym: bool) -> Array1<f64> {
        if m < 1 {
            return array![];
        } else if m % 2 == 1 {
            panic!("Kaiser-Bessel Derived windows are only defined for even number of points");
        }

        let mut kaiser_window =
            kaiser(((m as f64 / 2.0).floor() + 1.0) as u64, beta, true).to_vec();
        let mut last_sum = 0.0;
        for e in kaiser_window.iter_mut() {
            *e += last_sum;
            last_sum = *e;
        }
        let half_window = kaiser_window
            .iter()
            .enumerate()
            .take_while(|(i, _)| *i != kaiser_window.len() - 1)
            .map(|(_, a)| *a)
            .map(|a| (a / kaiser_window.last().unwrap()).sqrt());

        let mut w = half_window.collect::<Vec<_>>();
        let reversed: Vec<f64> = w.iter().copied().rev().collect();
        w.extend(reversed);

        w.into()
    }
    _kaiser_bessel_derived(m, beta, sym.into().unwrap_or(true))
}

pub fn gaussian(m: u64, std_dev: f64, sym: impl Into<Option<bool>>) -> Array1<f64> {
    fn _gaussian(m: u64, std_dev: f64, sym: bool) -> Array1<f64> {
        if_len_guard!(m);
        let (m, needs_trunc) = extend(m, sym);

        let n = Array1::range(0.0, m as _, 1.0) - (m as f64 - 1.0) / 2.0;
        let sig2 = 2.0 * std_dev * std_dev;

        let w = (-n.mapv(|a| a.powi(2)) / sig2).mapv(f64::exp);

        truncate(w, needs_trunc)
    }

    _gaussian(m, std_dev, sym.into().unwrap_or(true))
}

pub fn len_guards(m: u64) -> bool {
    m <= 1
}

pub fn extend(m: u64, sym: bool) -> (u64, bool) {
    if !sym {
        (m + 1, true)
    } else {
        (m, false)
    }
}

pub fn truncate(w: impl Into<Array1<f64>>, needs_trunc: bool) -> Array1<f64> {
    fn inner(w: Array1<f64>, needed: bool) -> Array1<f64> {
        if needed {
            let mut inner = w.to_vec();
            inner.remove(inner.len() - 1);
            inner.into()
        } else {
            w
        }
    }
    inner(w.into(), needs_trunc)
}

#[macro_export]
macro_rules! if_len_guard {
    ($m:ident) => {
        if len_guards($m) {
            return Array1::ones(($m as usize,));
        }
    };
}
