use std::ops::{Div, Mul};

use num::{Complex, Num, One, Zero};

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
    let z_prod: Complex<f64> = input
        .z
        .iter()
        .fold(Complex::one(), |acc, item| acc * (-item));
    let p_prod: Complex<f64> = input
        .p
        .iter()
        .fold(Complex::one(), |acc, item| acc * (-item));

    input.z = input.z.iter().map(|&a| wo / a).collect();
    input.p = input.p.iter().map(|&a| wo / a).collect();

    input.z.extend(vec![Complex::zero(); degree]);
    input.k = input.k * (z_prod / p_prod).re;
    println!("{}", input.k);
    input
}
pub fn lp2lp_zpk(mut input: Zpk, wo: f64) -> Zpk {
    let degree = relative_degree(&input);
    input.z = input.z.iter().map(|&a| a * wo).collect();
    input.p = input.p.iter().map(|&a| a * wo).collect();
    input.k = (input.k * wo).powi(degree as _);
    println!("{}", input.k);
    input
}
pub fn lp2bp_zpk(mut input: Zpk, wo: f64, bw: f64) -> Zpk {
    let degree = relative_degree(&input);
    input.z.iter_mut().for_each(|a| *a = *a * bw / 2.0);
    input.p.iter_mut().for_each(|a| *a = *a * bw / 2.0);

    let z_lp1 = input.z.iter().map(|&a| a + (a.powi(2) - wo.powi(2)).sqrt());
    let z_lp = z_lp1.chain(input.z.iter().map(|&a| a - (a.powi(2) - wo.powi(2)).sqrt()));
    input.z = z_lp.collect();

    let p_lp1 = input.p.iter().map(|&a| a + (a.powi(2) - wo.powi(2)).sqrt());
    let p_lp = p_lp1.chain(input.p.iter().map(|&a| a - (a.powi(2) - wo.powi(2)).sqrt()));
    input.p = p_lp.collect();

    input.z.extend(vec![Complex::zero(); degree]);

    input.k = input.k * bw.powi(degree as _);
    input
}
pub fn lp2bs_zpk(mut input: Zpk, wo: f64, bw: f64) -> Zpk {
    let degree = relative_degree(&input);

    let z_prod: Complex<f64> = input
        .z
        .iter()
        .fold(Complex::one(), |acc, item| acc * (-item));
    let p_prod: Complex<f64> = input
        .p
        .iter()
        .fold(Complex::one(), |acc, item| acc * (-item));

    input.z.iter_mut().for_each(|a| {
        *a = (bw / 2.0) / *a;
    });

    input.p.iter_mut().for_each(|a| {
        *a = (bw / 2.0) / *a;
    });

    let z_lp1 = input.z.iter().map(|&a| a + (a.powi(2) - wo.powi(2)).sqrt());
    let z_lp = z_lp1.chain(input.z.iter().map(|&a| a - (a.powi(2) - wo.powi(2)).sqrt()));
    input.z = z_lp.collect();

    let p_lp1 = input.p.iter().map(|&a| a + (a.powi(2) - wo.powi(2)).sqrt());
    let p_lp = p_lp1.chain(input.p.iter().map(|&a| a - (a.powi(2) - wo.powi(2)).sqrt()));
    input.p = p_lp.collect();

    input.z.extend(vec![Complex::new(0.0, 1.0) * wo; degree]);
    input.z.extend(vec![Complex::new(0.0, -1.0) * wo; degree]);

    input.k = input.k * (z_prod / p_prod).re;

    input
}
