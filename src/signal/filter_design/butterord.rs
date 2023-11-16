use ndarray::{array, Array1};
use num::{Float, NumCast};
use thiserror::Error;

use crate::signal::{band_filter::GenericOrdBandFilter, GenericSampling};

use super::{OrdCompute, OrdResult};

pub struct ButterOrd<T> {
    pub band_filter: GenericOrdBandFilter<T>,
    pub gpass: T,
    pub gstop: T,
    pub sampling: GenericSampling<T>,
}

impl<T: Float> OrdCompute<T> for ButterOrd<T> {
    fn compute_order(&self) -> Result<OrdResult<T>, crate::signal::error::Error> {
        Ok(
            buttord(self.band_filter, self.gpass, self.gstop, self.sampling)
                .map_err(super::Error::from)?,
        )
    }
}

pub fn buttord<T: Float>(
    band_filter: GenericOrdBandFilter<T>,
    gpass: T,
    gstop: T,
    sampling: GenericSampling<T>,
) -> Result<OrdResult<T>, Error> {
    _validate_gpass_gstop(gpass, gstop)?;
    let band_filter = _validate_wp_ws(band_filter, sampling);
    let band_filter = _pre_warp(band_filter, sampling);
    let (nat, band_filter) = _find_nat_freq(band_filter, gpass, gstop);

    let g_stop = T::from(10.0)
        .unwrap()
        .powf(T::from(0.1).unwrap() * gstop.abs());
    let g_pass = T::from(10.0)
        .unwrap()
        .powf(T::from(0.1).unwrap() * gpass.abs());

    // int(ceil(log10((GSTOP - 1.0) / (GPASS - 1.0)) / (2 * log10(nat))))
    let ord = <f64 as NumCast>::from(
        (((g_stop - T::one()) / (g_pass - T::one())).log10() / (T::from(2).unwrap() * nat.log10()))
            .ceil(),
    )
    .unwrap();

    let w0 = (g_pass - T::one()).powf(T::from(-1.0 / (2.0 * ord)).unwrap());
    let ord = ord as u32;

    use crate::signal::band_filter::GenericOrdBandFilterType as S;
    let band_filter = match band_filter.as_ref() {
        S::Lowpass { .. } => band_filter * w0,
        S::Highpass { .. } => band_filter / w0,
        S::Bandpass {
            wp_low,
            wp_high,
            ws_low,
            ws_high,
        } => todo!(),
        S::Bandstop {
            wp_low,
            wp_high,
            ws_low,
            ws_high,
        } => todo!(),
    };
    todo!()
}

fn _validate_gpass_gstop<T: Float>(gpass: T, gstop: T) -> Result<(), Error> {
    if gpass <= T::zero() {
        return Err(Error::BadGPass(gpass.to_f64().unwrap()));
    }
    if gstop <= T::zero() {
        return Err(Error::BadGStop(gstop.to_f64().unwrap()));
    }
    if gpass > gstop {
        return Err(Error::BadGpassAndGstop {
            gstop: gstop.to_f64().unwrap(),
            gpass: gpass.to_f64().unwrap(),
        });
    }
    Ok(())
}

fn _validate_wp_ws<T: Float>(
    mut band_filter: GenericOrdBandFilter<T>,
    sampling: GenericSampling<T>,
) -> GenericOrdBandFilter<T> {
    if let GenericSampling::Digital { fs } = &sampling {
        band_filter = (band_filter * T::from(2.0).unwrap()) / *fs;
    }

    band_filter
}

fn _pre_warp<T: Float>(
    mut band_filter: GenericOrdBandFilter<T>,
    sampling: GenericSampling<T>,
) -> GenericOrdBandFilter<T> {
    use std::f64::consts::PI;

    if !sampling.is_analog() {
        band_filter = ((band_filter * T::from(PI).unwrap()) / T::from(2.0).unwrap()).tan();
    }

    band_filter
}

fn _find_nat_freq<T: Float>(
    band_filter: GenericOrdBandFilter<T>,
    gpass: T,
    gstop: T,
) -> (T, GenericOrdBandFilter<T>) {
    let nat = match *band_filter.as_ref() {
        crate::signal::band_filter::GenericOrdBandFilterType::Lowpass { wp, ws } => wp / ws,
        crate::signal::band_filter::GenericOrdBandFilterType::Highpass { wp, ws } => ws / wp,
        crate::signal::band_filter::GenericOrdBandFilterType::Bandpass {
            wp_low,
            wp_high,
            ws_low,
            ws_high,
        } => {
            let nat_1 = (ws_low.powi(2) - wp_low * wp_high) / (ws_low * (wp_low - wp_high));
            let nat_2 = (ws_high.powi(2) - wp_low * wp_high) / (ws_high * (wp_low - wp_high));
            nat_1.abs().min(nat_2.abs())
        }
        crate::signal::band_filter::GenericOrdBandFilterType::Bandstop {
            wp_low,
            wp_high,
            ws_low,
            ws_high,
        } => unimplemented!(),
    };
    (nat, band_filter)
}

#[derive(Debug, Error)]
pub enum Error {
    #[error("gpass should be larger than 0.0, received {0}")]
    BadGPass(f64),
    #[error("gstop should be larger than 0.0, received {0}")]
    BadGStop(f64),
    #[error("gpass should be smaller than gstop, received: gpass {gpass}, gstop {gstop}")]
    BadGpassAndGstop { gpass: f64, gstop: f64 },
}
