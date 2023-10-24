use crate::signal::band_filter::GenericBandFilter;
use crate::signal::output_type::{Ba, FilterOutput, GenericBa, GenericFilterOutput};
use crate::signal::{BandFilter, GenericSampling, Sampling};
use crate::special::sinc;
use ndarray::{array, Array1};
use num::Float;

use super::windows::{get_window, WindowType};
use super::{kaiser_atten, kaiser_beta, GenericFIRFilterSettings};

pub struct Firwin1Filter<T> {
    pub settings: GenericFIRFilterSettings<T>,
}

impl Firwin1Filter<f64> {
    pub fn firwin(self) -> FilterOutput {
        let GenericFIRFilterSettings {
            numtaps,
            cutoff,
            width,
            window,
            scale,
            sampling,
        } = self.settings;

        firwin(numtaps, cutoff, width, window, scale, sampling)
    }
}

pub fn firwin<T: Float>(
    numtaps: i64,
    cutoff: GenericBandFilter<T>,
    width: Option<T>,
    mut window: WindowType<T>,
    scale: bool,
    sampling: GenericSampling<T>,
) -> GenericFilterOutput<T> {
    let nyq = match sampling {
        GenericSampling::Digital { fs } => T::from(0.5).unwrap() * fs,
        GenericSampling::Analog => T::one(),
    };
    let cutoff = cutoff / nyq;

    let pass_zero = cutoff.pass_zero();

    let pass_nyquist = cutoff.pass_nyquist(pass_zero);

    if pass_nyquist && numtaps % 2 == 0 {
        panic!("A filter with an even number of coefficients must have zero response at the Nyquist frequency.");
    }

    if let Some(width) = width {
        let atten = kaiser_atten(numtaps, width);
        let beta = kaiser_beta(atten);
        window = WindowType::Kaiser { beta };
    }

    let mut cutoff = cutoff.to_vec();
    if pass_zero {
        cutoff.insert(0, T::zero());
    }
    if pass_nyquist {
        cutoff.push(T::one());
    }

    let size = cutoff.len();
    if size % 2 != 0 {
        panic!("");
    }

    let cutoff = Array1::from_vec(cutoff);
    let bands = cutoff.into_shape((size / 2, 2)).unwrap();

    let alpha = 0.5 * ((numtaps as f64) - 1.0);

    let m = (0..numtaps)
        .map(|a| T::from((a as f64) - alpha).unwrap())
        .collect::<Array1<_>>();

    let mut h = Array1::from_vec(vec![T::zero(); numtaps as usize]);

    for row in bands.rows() {
        let left = row[0];
        let right = row[1];

        h = h + &(sinc(m.mapv(|a| a * right)).mapv(|a| a * right));
        h = h - &(sinc(m.mapv(|a| a * left)).mapv(|a| a * left));
    }
    let win = get_window(window, numtaps as _, false);
    h = h * win;

    if scale {
        let first = bands.rows().into_iter().next().unwrap();
        let (left, right) = (first[0], first[1]);

        let scale_frequency = if left == T::zero() {
            T::zero()
        } else if right == T::one() {
            T::one()
        } else {
            T::from(0.5).unwrap() * (left + right)
        };

        let c = (m.mapv(|a| a * T::from(std::f64::consts::PI).unwrap() * scale_frequency))
            .mapv(|a| a.cos());
        let s = (c * &h).sum();
        h = h.mapv(|a| a / s);
    }

    GenericFilterOutput::Ba(GenericBa {
        a: array![T::one()].mapv(Into::into),
        b: h.mapv(Into::into),
    })
}
