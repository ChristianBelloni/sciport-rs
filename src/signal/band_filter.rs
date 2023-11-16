use super::{output_type::GenericZpk, tools::relative_degree};
use ndarray::{array, concatenate, Array1, ArrayView, Axis};
use num::{Complex, Float, Num, NumCast, Zero};
use std::{
    borrow::Cow,
    ops::{Add, Div, Mul, Sub},
};
use thiserror::Error;

pub type BandFilter = GenericBandFilter<f64>;

#[derive(Debug, Clone, Copy)]
pub enum GenericBandFilter<T> {
    Highpass(T),
    Lowpass(T),
    Bandpass { low: T, high: T },
    Bandstop { low: T, high: T },
}

impl<T> GenericBandFilter<T> {
    pub fn cast<K>(self) -> GenericBandFilter<K>
    where
        K: From<T>,
    {
        match self {
            GenericBandFilter::Highpass(data) => GenericBandFilter::Highpass(data.into()),
            GenericBandFilter::Lowpass(data) => GenericBandFilter::Lowpass(data.into()),
            GenericBandFilter::Bandpass { low, high } => GenericBandFilter::Bandpass {
                low: low.into(),
                high: high.into(),
            },
            GenericBandFilter::Bandstop { low, high } => GenericBandFilter::Bandstop {
                low: low.into(),
                high: high.into(),
            },
        }
    }

    pub fn cast_with_fn<K>(self, f: impl Fn(T) -> K) -> GenericBandFilter<K> {
        match self {
            GenericBandFilter::Highpass(data) => GenericBandFilter::Highpass(f(data)),
            GenericBandFilter::Lowpass(data) => GenericBandFilter::Lowpass(f(data)),
            GenericBandFilter::Bandpass { low, high } => GenericBandFilter::Bandpass {
                low: f(low),
                high: f(high),
            },
            GenericBandFilter::Bandstop { low, high } => GenericBandFilter::Bandstop {
                low: f(low),
                high: f(high),
            },
        }
    }
}

impl<T> GenericBandFilter<T>
where
    T: Float,
{
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
            Self::Highpass(_) | Self::Lowpass(_) => 1,
            Self::Bandstop { low: _, high: _ } | Self::Bandpass { low: _, high: _ } => 2,
        }
    }
    pub fn to_vec(self) -> Vec<T> {
        match self {
            Self::Highpass(f) | Self::Lowpass(f) => vec![f],
            Self::Bandstop { low, high } | Self::Bandpass { low, high } => {
                vec![low, high]
            }
        }
    }
    pub fn to_array(self) -> Array1<T> {
        match self {
            Self::Highpass(f) | Self::Lowpass(f) => array![f],
            Self::Bandstop { low, high } | Self::Bandpass { low, high } => {
                array![low, high]
            }
        }
    }

    pub fn pass_zero(&self) -> bool {
        match self {
            Self::Bandpass { low: _, high: _ } | Self::Highpass(_) => false,
            Self::Bandstop { low: _, high: _ } | Self::Lowpass(_) => true,
        }
    }

    pub fn pass_nyquist(&self, pass_zero: bool) -> bool {
        ((self.size() & 1) == 1) ^ pass_zero
    }
}

impl<T: Num + Copy, R: Num + Copy> Mul<R> for GenericBandFilter<T>
where
    T: Mul<R, Output = T>,
{
    type Output = Self;
    fn mul(self, rhs: R) -> Self::Output {
        match self {
            Self::Highpass(data) => Self::Highpass(data * rhs),
            Self::Lowpass(data) => Self::Lowpass(data * rhs),
            Self::Bandpass { low, high } => Self::Bandpass {
                low: low * rhs,
                high: high * rhs,
            },
            Self::Bandstop { low, high } => Self::Bandstop {
                low: low * rhs,
                high: high * rhs,
            },
        }
    }
}

impl<T: Num + Copy> Div<T> for GenericBandFilter<T> {
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        match self {
            Self::Highpass(data) => Self::Highpass(data / rhs),
            Self::Lowpass(data) => Self::Lowpass(data / rhs),
            Self::Bandpass { low, high } => Self::Bandpass {
                low: low / rhs,
                high: high / rhs,
            },
            Self::Bandstop { low, high } => Self::Bandstop {
                low: low / rhs,
                high: high / rhs,
            },
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GenericOrdBandFilter<T>(GenericOrdBandFilterType<T>);

impl<T: std::ops::Deref> std::ops::Deref for GenericOrdBandFilter<T> {
    type Target = GenericOrdBandFilterType<T>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> AsRef<GenericOrdBandFilterType<T>> for GenericOrdBandFilter<T> {
    fn as_ref(&self) -> &GenericOrdBandFilterType<T> {
        &self.0
    }
}

#[derive(Debug, Error)]
pub enum Error {
    #[error("{0}")]
    Validation(Cow<'static, str>),
}

impl<T: Float> GenericOrdBandFilter<T> {
    pub fn lowpass(wp: T, ws: T) -> Result<Self, Error> {
        if wp > ws {
            Err(Error::Validation(Cow::from(
                "for lowpass filter wp must be smaller than ws",
            )))?;
        }
        Ok(Self(GenericOrdBandFilterType::Lowpass { wp, ws }))
    }

    pub fn highpass(wp: T, ws: T) -> Result<Self, Error> {
        if wp < ws {
            Err(Error::Validation(Cow::from(
                "for highpass filter ws must be smaller than wp",
            )))?;
        }
        Ok(Self(GenericOrdBandFilterType::Lowpass { wp, ws }))
    }

    pub fn bandpass(wp_low: T, wp_high: T, ws_low: T, ws_high: T) -> Result<Self, Error> {
        if wp_low < ws_low {
            Err(Error::Validation(Cow::from(
                "for bandpass filter ws_low must be smaller than wp_low",
            )))?;
        }
        Ok(Self(GenericOrdBandFilterType::Bandpass {
            wp_low,
            wp_high,
            ws_low,
            ws_high,
        }))
    }

    pub fn bandstop(wp_low: T, wp_high: T, ws_low: T, ws_high: T) -> Result<Self, Error> {
        if ws_low < wp_low {
            Err(Error::Validation(Cow::from(
                "for bandstop filter wp_low must be smaller than ws_low",
            )))?;
        }
        Ok(Self(GenericOrdBandFilterType::Bandpass {
            wp_low,
            wp_high,
            ws_low,
            ws_high,
        }))
    }

    pub fn tan(self) -> Self {
        use GenericOrdBandFilterType as S;
        let inner = match self.0 {
            S::Lowpass { wp, ws } => S::Lowpass {
                wp: wp.tan(),
                ws: ws.tan(),
            },
            S::Highpass { wp, ws } => S::Highpass {
                wp: wp.tan(),
                ws: ws.tan(),
            },
            S::Bandpass {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandpass {
                wp_low: wp_low.tan(),
                wp_high: wp_high.tan(),
                ws_low: ws_low.tan(),
                ws_high: ws_high.tan(),
            },
            S::Bandstop {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandstop {
                wp_low: wp_low.tan(),
                wp_high: wp_high.tan(),
                ws_low: ws_low.tan(),
                ws_high: ws_high.tan(),
            },
        };
        Self(inner)
    }
}

#[derive(Debug, Clone, Copy)]
pub enum GenericOrdBandFilterType<T> {
    Lowpass {
        wp: T,
        ws: T,
    },
    Highpass {
        wp: T,
        ws: T,
    },
    Bandpass {
        wp_low: T,
        wp_high: T,
        ws_low: T,
        ws_high: T,
    },
    Bandstop {
        wp_low: T,
        wp_high: T,
        ws_low: T,
        ws_high: T,
    },
}

impl<T: Float> Div<T> for GenericOrdBandFilter<T> {
    type Output = Self;

    fn div(self, rhs: T) -> Self::Output {
        use GenericOrdBandFilterType as S;

        let inner = match self.0 {
            S::Lowpass { wp, ws } => S::Lowpass {
                wp: wp / rhs,
                ws: ws / rhs,
            },
            S::Highpass { wp, ws } => S::Highpass {
                wp: wp / rhs,
                ws: ws / rhs,
            },
            S::Bandpass {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandpass {
                wp_low: wp_low / rhs,
                wp_high: wp_high / rhs,
                ws_low: ws_low / rhs,
                ws_high: ws_high / rhs,
            },
            S::Bandstop {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandstop {
                wp_low: wp_low / rhs,
                wp_high: wp_high / rhs,
                ws_low: ws_low / rhs,
                ws_high: ws_high / rhs,
            },
        };
        Self(inner)
    }
}
impl<T: Float> Sub<T> for GenericOrdBandFilter<T> {
    type Output = Self;

    fn sub(self, rhs: T) -> Self::Output {
        use GenericOrdBandFilterType as S;

        let inner = match self.0 {
            S::Lowpass { wp, ws } => S::Lowpass {
                wp: wp - rhs,
                ws: ws - rhs,
            },
            S::Highpass { wp, ws } => S::Highpass {
                wp: wp - rhs,
                ws: ws - rhs,
            },
            S::Bandpass {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandpass {
                wp_low: wp_low - rhs,
                wp_high: wp_high - rhs,
                ws_low: ws_low - rhs,
                ws_high: ws_high - rhs,
            },
            S::Bandstop {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandstop {
                wp_low: wp_low - rhs,
                wp_high: wp_high - rhs,
                ws_low: ws_low - rhs,
                ws_high: ws_high - rhs,
            },
        };
        Self(inner)
    }
}
impl<T: Float> Add<T> for GenericOrdBandFilter<T> {
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        use GenericOrdBandFilterType as S;

        let inner = match self.0 {
            S::Lowpass { wp, ws } => S::Lowpass {
                wp: wp + rhs,
                ws: ws + rhs,
            },
            S::Highpass { wp, ws } => S::Highpass {
                wp: wp + rhs,
                ws: ws + rhs,
            },
            S::Bandpass {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandpass {
                wp_low: wp_low + rhs,
                wp_high: wp_high + rhs,
                ws_low: ws_low + rhs,
                ws_high: ws_high + rhs,
            },
            S::Bandstop {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandstop {
                wp_low: wp_low + rhs,
                wp_high: wp_high + rhs,
                ws_low: ws_low + rhs,
                ws_high: ws_high + rhs,
            },
        };
        Self(inner)
    }
}

impl<T: Float> Mul<T> for GenericOrdBandFilter<T> {
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        use GenericOrdBandFilterType as S;

        let inner = match self.0 {
            S::Lowpass { wp, ws } => S::Lowpass {
                wp: wp * rhs,
                ws: ws * rhs,
            },
            S::Highpass { wp, ws } => S::Highpass {
                wp: wp * rhs,
                ws: ws * rhs,
            },
            S::Bandpass {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandpass {
                wp_low: wp_low * rhs,
                wp_high: wp_high * rhs,
                ws_low: ws_low * rhs,
                ws_high: ws_high * rhs,
            },
            S::Bandstop {
                wp_low,
                wp_high,
                ws_low,
                ws_high,
            } => S::Bandstop {
                wp_low: wp_low * rhs,
                wp_high: wp_high * rhs,
                ws_low: ws_low * rhs,
                ws_high: ws_high * rhs,
            },
        };
        Self(inner)
    }
}

pub fn lp2bf_zpk<T>(input: GenericZpk<T>, wo: GenericBandFilter<T>) -> GenericZpk<T>
where
    T: Float + Clone,
{
    match wo {
        GenericBandFilter::Lowpass(wo) => lp2lp_zpk(input, wo),
        GenericBandFilter::Highpass(wo) => lp2hp_zpk(input, wo),
        GenericBandFilter::Bandpass { low, high } => {
            let bw = high - low;
            let wo = (low * high).sqrt();
            lp2bp_zpk(input, wo, bw)
        }
        GenericBandFilter::Bandstop { low, high } => {
            let bw = high - low;
            let wo = (low * high).sqrt();
            lp2bs_zpk(input, wo, bw)
        }
    }
}

pub fn lp2hp_zpk<T>(mut input: GenericZpk<T>, wo: T) -> GenericZpk<T>
where
    T: Float + Clone,
{
    let degree = relative_degree(&input);

    let big_z = &input.z;
    let big_p = &input.p;

    let z_prod = big_z.map(|a| -a).product();
    let p_prod = big_p.map(|a| -a).product();

    input
        .z
        .mapv_inplace(|a| <num::Complex<T> as From<T>>::from(wo) / a);
    input
        .p
        .mapv_inplace(|a| <num::Complex<T> as From<T>>::from(wo) / a);

    let zeros = vec![Complex::zero(); degree];
    input
        .z
        .append(Axis(0), ArrayView::from(zeros.as_slice()))
        .unwrap();

    let factor = (z_prod / p_prod).re;

    input.k = input.k * factor;
    if input.k.is_nan() {
        println!("lp2hp nan")
    }
    input
}

pub fn lp2lp_zpk<T>(mut input: GenericZpk<T>, wo: T) -> GenericZpk<T>
where
    T: Float,
{
    let degree = relative_degree(&input);
    input.z.mapv_inplace(|a| a * wo);
    input.p.mapv_inplace(|a| a * wo);
    let res = input.k * wo.powi(degree as _);
    res.is_nan();
    input.k = res;
    input
}

pub fn lp2bp_zpk<T>(input: GenericZpk<T>, wo: T, bw: T) -> GenericZpk<T>
where
    T: Float,
{
    let degree = relative_degree(&input);
    let GenericZpk { z, p, k } = input;
    let two = T::one() + T::one();
    let z_lp = z.mapv(|a| a * bw / two);
    let p_lp = p.mapv(|a| a * bw / two);

    let z_hp_left = &z_lp + (z_lp.mapv(|a| a.powf(two) - wo.powf(two))).mapv(Complex::sqrt);
    let z_hp_right = &z_lp - (z_lp.mapv(|a| a.powf(two) - wo.powf(two))).mapv(Complex::sqrt);

    let p_hp_left = &p_lp + (p_lp.mapv(|a| a.powf(two) - wo.powf(two))).mapv(Complex::sqrt);
    let p_hp_right = &p_lp - (p_lp.mapv(|a| a.powf(two) - wo.powf(two))).mapv(Complex::sqrt);

    let z_bp = concatenate![Axis(0), z_hp_left, z_hp_right];
    let p_bp = concatenate![Axis(0), p_hp_left, p_hp_right];

    let z_bp = concatenate![Axis(0), z_bp, Array1::zeros(degree)];

    let k_bp = k * bw.powi(degree as _);

    let p_bp = p_bp.mapv(|mut a| {
        if a.im == T::neg_zero() {
            a.im = T::zero()
        }
        a
    });

    GenericZpk {
        z: z_bp,
        p: p_bp,
        k: k_bp,
    }
}

pub fn lp2bs_zpk<T>(input: GenericZpk<T>, wo: T, bw: T) -> GenericZpk<T>
where
    T: Float,
{
    let degree = relative_degree(&input);
    let GenericZpk { z, p, k } = input;

    let bw_half = bw / T::from(2).unwrap();
    let bw_half = Complex::new(bw_half, T::zero());

    let z_hp = z.mapv(|a| bw_half / a);
    let p_hp = p.mapv(|a| bw_half / a);

    let z_bs_left = &z_hp + (z_hp.mapv(|a| a.powi(2) - wo.powi(2))).mapv(Complex::sqrt);
    let z_bs_right = &z_hp - (z_hp.mapv(|a| a.powi(2) - wo.powi(2))).mapv(Complex::sqrt);

    let p_bs_left = &p_hp + (p_hp.mapv(|a| a.powi(2) - wo.powi(2))).mapv(Complex::sqrt);
    let p_bs_right = &p_hp - (p_hp.mapv(|a| a.powi(2) - wo.powi(2))).mapv(Complex::sqrt);

    let dbg_zbs_left = z_bs_left.mapv(|a| {
        Complex::new(
            <f64 as NumCast>::from(a.re).unwrap(),
            <f64 as NumCast>::from(a.im).unwrap(),
        )
    });

    let dbg_zbs_right = z_bs_right.mapv(|a| {
        Complex::new(
            <f64 as NumCast>::from(a.re).unwrap(),
            <f64 as NumCast>::from(a.im).unwrap(),
        )
    });

    println!("dbg z {dbg_zbs_left:?}");
    println!("dbg z {dbg_zbs_right:?}");

    let z_bs = concatenate![Axis(0), z_bs_left, z_bs_right];
    let p_bs = concatenate![Axis(0), p_bs_left, p_bs_right];

    let dbg_zbs = z_bs.mapv(|a| {
        Complex::new(
            <f64 as NumCast>::from(a.re).unwrap(),
            <f64 as NumCast>::from(a.im).unwrap(),
        )
    });
    println!("dbg z {dbg_zbs:?}");
    let z_bs = concatenate![Axis(0), z_bs, Array1::from_elem(degree, Complex::i() * wo)];
    let z_bs = concatenate![
        Axis(0),
        z_bs,
        Array1::from_elem(degree, Complex::new(T::zero(), -T::one()) * wo)
    ];

    println!("degree {degree} p length: {}", p.len());
    let factor = match z.len().cmp(&p.len()) {
        std::cmp::Ordering::Less => {
            let t_z = concatenate![Axis(0), -&z, Array1::ones(p.len() - z.len())];
            (t_z / -&p).product()
        }
        std::cmp::Ordering::Equal => (-&z / -&p).product(),
        std::cmp::Ordering::Greater => {
            let t_p = concatenate![Axis(0), -&p, Array1::ones(z.len() - p.len())];
            (&-z / t_p).product()
        }
    };

    let k_bs = k * factor.re;

    GenericZpk {
        z: z_bs,
        p: p_bs,
        k: k_bs,
    }
}
