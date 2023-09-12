use std::ops::{Div, Mul};

use ndarray::{array, concatenate, Array, Array1, ArrayView, Axis};
use num::{traits::Pow, Complex, Num, Zero};

use super::{output_type::Zpk, tools::relative_degree};

#[derive(Debug, Clone, Copy)]
pub enum BandFilter {
    Highpass(f64),
    Lowpass(f64),
    Bandpass { low: f64, high: f64 },
    Bandstop { low: f64, high: f64 },
}

impl BandFilter {
    pub fn tan(self) -> Self {
        match self {
            Self::Highpass(data) => Self::Highpass(data.tan()),
            Self::Lowpass(data) => Self::Lowpass(data.tan()),
            Self::Bandpass { low, high } => Self::Bandpass {
                low: low.tan(),
                high: high.tan(),
            },
            Self::Bandstop { low, high } => Self::Bandstop {
                low: low.tan(),
                high: high.tan(),
            },
        }
    }

    pub fn size(&self) -> u8 {
        match self {
            BandFilter::Highpass(_) | BandFilter::Lowpass(_) => 1,
            BandFilter::Bandstop { low: _, high: _ } | BandFilter::Bandpass { low: _, high: _ } => {
                2
            }
        }
    }
    pub fn to_vec(self) -> Vec<f64> {
        match self {
            BandFilter::Highpass(f) | BandFilter::Lowpass(f) => vec![f],
            BandFilter::Bandstop { low, high } | BandFilter::Bandpass { low, high } => {
                vec![low, high]
            }
        }
    }
    pub fn to_array(self) -> Array1<f64> {
        match self {
            BandFilter::Highpass(f) | BandFilter::Lowpass(f) => array![f],
            BandFilter::Bandstop { low, high } | BandFilter::Bandpass { low, high } => {
                array![low, high]
            }
        }
    }
}

impl<T: Num + Copy> Mul<T> for BandFilter
where
    f64: From<T>,
    T: From<f64>,
{
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        match self {
            Self::Highpass(data) => Self::Highpass((rhs.mul(data.into())).into()),
            Self::Lowpass(data) => Self::Lowpass((rhs.mul(data.into())).into()),
            Self::Bandpass { low, high } => Self::Bandpass {
                low: (rhs.mul(low.into())).into(),
                high: (rhs.mul(high.into())).into(),
            },
            Self::Bandstop { low, high } => Self::Bandstop {
                low: (rhs.mul(low.into())).into(),
                high: (rhs.mul(high.into())).into(),
            },
        }
    }
}
impl<T: Num + Copy> Div<T> for BandFilter
where
    f64: From<T>,
    T: From<f64>,
{
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        match self {
            Self::Highpass(data) => Self::Highpass(data / Into::<f64>::into(rhs)),
            Self::Lowpass(data) => Self::Lowpass(data / Into::<f64>::into(rhs)),
            Self::Bandpass { low, high } => Self::Bandpass {
                low: low / Into::<f64>::into(rhs),
                high: high / Into::<f64>::into(rhs),
            },
            Self::Bandstop { low, high } => Self::Bandstop {
                low: low / Into::<f64>::into(rhs),
                high: high / Into::<f64>::into(rhs),
            },
        }
    }
}

pub fn lp2bf_zpk(input: Zpk, wo: BandFilter) -> Zpk {
    match wo {
        BandFilter::Lowpass(wo) => lp2lp_zpk(input, wo),
        BandFilter::Highpass(wo) => lp2hp_zpk(input, wo),
        BandFilter::Bandpass { low, high } => {
            let bw = high - low;
            let wo = (high * low).sqrt();
            lp2bp_zpk(input, wo, bw)
        }
        BandFilter::Bandstop { low, high } => {
            let bw = high - low;
            let wo = (high * low).sqrt();
            lp2bs_zpk(input, wo, bw)
        }
    }
}

pub fn lp2hp_zpk(mut input: Zpk, wo: f64) -> Zpk {
    let degree = relative_degree(&input);

    let z_prod = input.z.map(|a| -a).product();
    let p_prod = input.p.map(|a| -a).product();

    input.z.mapv_inplace(|a| wo / a);
    input.p.mapv_inplace(|a| wo / a);

    let zeros = vec![Complex::zero(); degree];
    input
        .z
        .append(Axis(0), ArrayView::from(zeros.as_slice()))
        .unwrap();

    input.k *= (z_prod / p_prod).re;
    input
}

pub fn lp2lp_zpk(mut input: Zpk, wo: f64) -> Zpk {
    let degree = relative_degree(&input);
    input.z.mapv_inplace(|a| a * wo);
    input.p.mapv_inplace(|a| a * wo);
    input.k = input.k * wo.powi(degree as _);
    input
}
pub fn lp2bp_zpk(mut input: Zpk, wo: f64, bw: f64) -> Zpk {
    let degree = relative_degree(&input);
    let wo = Complex::from(wo);
    let bw = Complex::from(bw);

    let z_lp = input.z.mapv(|a| a * Complex::from(bw) / Complex::from(2.0));
    let p_lp = input.p.mapv(|a| a * Complex::from(bw) / Complex::from(2.0));

    let z_bp1 = z_lp.mapv(|a| a + Complex::from(a.pow(2.0) - wo.pow(2.0)).sqrt());

    let z_bp2 = z_lp.mapv(|a| a - Complex::from(a.pow(2.0) - wo.pow(2.0)).sqrt());

    let z_bp = concatenate![Axis(0), z_bp1, z_bp2];

    let p_bp1 = p_lp.mapv(|a| a + Complex::sqrt(Complex::from(a.powf(2.0) - wo.powf(2.0))));

    let p_bp2 = p_lp.mapv(|a| a - Complex::sqrt(Complex::from(a.powf(2.0) - wo.powf(2.0))));

    let p_bp = concatenate![Axis(0), p_bp1, p_bp2];
    let z_bp = concatenate![
        Axis(0),
        z_bp,
        Array::from_shape_vec(degree, vec![Complex::<f64>::zero(); degree]).unwrap()
    ];

    input.k *= bw.powi(degree as _).norm();
    input.z = z_bp;
    input.p = p_bp;
    input
}

pub fn lp2bs_zpk(mut input: Zpk, wo: f64, bw: f64) -> Zpk {
    let degree = relative_degree(&input);

    let bw_half = bw / 2.0;
    let z_hp = input.z.mapv(|a: Complex<f64>| bw_half / a).to_owned();
    let mut p_hp = Complex::from(bw_half) / input.p.clone();

    // this is weird, i have no idea why but removing this breaks the conversion to BS
    // Lord forgive me
    p_hp.mapv_inplace(|mut a| {
        if a.im == 0.0 {
            a.im = -0.0;
        }
        if a.re == 0.0 {
            a.re = -0.0;
        }
        a
    });

    let z_bs1 = z_hp.map(|a| a + (a.powi(2) - wo.powi(2)).sqrt());

    let z_bs2 = z_hp.map(|a| a - (a.powi(2) - wo.powi(2)).sqrt());

    let z_bs = concatenate![Axis(0), z_bs1, z_bs2];

    let wo_squared = Complex::from(wo.pow(2));

    let p_bs1 =
        p_hp.clone() + (p_hp.mapv(|a| Complex::powi(&a, 2)) - wo_squared).mapv(Complex::sqrt);
    let p_bs2 =
        p_hp.clone() - (p_hp.mapv(|a| Complex::powi(&a, 2)) - wo_squared).mapv(Complex::sqrt);

    let p_bs = concatenate![Axis(0), p_bs1, p_bs2];

    let z_bs = concatenate![
        Axis(0),
        z_bs,
        Array::from_shape_vec(degree, vec![Complex::<f64>::new(0.0, wo); degree]).unwrap()
    ];
    let z_bs = concatenate![
        Axis(0),
        z_bs,
        Array::from_shape_vec(degree, vec![Complex::<f64>::new(0.0, -wo); degree]).unwrap()
    ];

    let z_prod = input.z.mapv(|a| -a).product();
    let p_prod = input.p.map(|a| -a).product();
    let factor = (z_prod / p_prod).re;
    input.k *= factor;
    input.z = z_bs;
    input.p = p_bs;
    input
}
