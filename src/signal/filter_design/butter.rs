use std::marker::PhantomData;

use ndarray::{array, Array1};
use num::{complex::Complex64, Complex};

use crate::signal::{
    band_filter::BandFilter,
    iir_filter,
    output_type::{Ba, DesiredFilterOutput, FilterOutput, Zpk},
    Analog,
};

/// Butterworth digital and analog filter design.
///
/// Design an Nth-order digital or analog Butterworth filter and return the filter coefficients.
///
/// # Notes
///
/// The Butterworth filter has maximally flat frequency response in the passband.
///
/// See [`buttap`] for implementation details and references.
///
pub struct ButterFilter<T> {
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    cache: Option<T>,
}

impl<T> ButterFilter<T> {
    pub fn new(order: u32, band_filter: BandFilter, analog: Analog) -> Self {
        Self {
            order,
            band_filter,
            analog,
            cache: None,
        }
    }
}

crate::impl_iir!(
    ButterFilter<Zpk>,
    Zpk,
    self,
    ButterFilterStandalone::<Zpk>::filter(self.order, self.band_filter, self.analog)
);

crate::impl_iir!(
    ButterFilter<Ba>,
    Ba,
    self,
    ButterFilterStandalone::<Ba>::filter(self.order, self.band_filter, self.analog)
);

/// Butterworth digital and analog filter design.
///
/// Design an Nth-order digital or analog Butterworth filter and return the filter coefficients.
///
/// # Notes
///
/// The Butterworth filter has maximally flat frequency response in the passband.
///
/// See [`buttap`] for implementation details and references.
///
pub struct ButterFilterStandalone<T>(PhantomData<T>);

impl ButterFilterStandalone<Zpk> {
    /// Butterworth digital and analog filter design.
    ///
    /// Design an Nth-order digital or analog Butterworth filter and return the filter coefficients.
    ///
    /// # Parameters
    ///  - order : Order of the filter
    ///  - band: [`BandFilter`] to apply, for analog filters band is expressed as an angular</br>
    /// frequency (rad/s).
    /// for digital filters band is in the same units as [`Analog`]
    ///  - analog: Analog or Digital filter selection, when digital is selected<br/>
    ///  a sampling rate is required
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog) -> Zpk {
        let filter = butter_filter(order, band_filter, analog, DesiredFilterOutput::Zpk);
        filter.zpk()
    }
}

impl ButterFilterStandalone<Ba> {
    /// Butterworth digital and analog filter design.
    ///
    /// Design an Nth-order digital or analog Butterworth filter and return the filter coefficients.
    ///
    /// # Parameters
    ///  - order : Order of the filter
    ///  - band: [`BandFilter`] to apply, for analog filters band is expressed as an angular</br>
    /// frequency (rad/s).
    /// for digital filters band is in the same units as [`Analog`]
    ///  - analog: Analog or Digital filter selection, when digital is selected<br/>
    ///  a sampling rate is required
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog) -> Ba {
        let filter = butter_filter(order, band_filter, analog, DesiredFilterOutput::Ba);
        filter.ba()
    }
}

pub(crate) fn butter_filter(
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let proto = buttap(order);
    iir_filter(proto, order, band_filter, analog, desired_output)
}

pub fn buttap(order: u32) -> Zpk {
    let order = order as i32;
    let z = array![];
    let mut range: Array1<Complex64> = ((-order + 1)..order)
        .step_by(2)
        .map(|a| (a as f64).into())
        .collect();

    fn make_iteration(item: Complex64, order: i32) -> Complex<f64> {
        let numerator = Complex::new(0.0, 1.0) * std::f64::consts::PI * item;
        let denominator = 2.0 * order as f64;

        let mut temp = numerator / denominator;
        if temp.im == 0.0 || temp.im == -0.0 {
            temp.im = 0.0;
        }
        -temp.exp()
    }

    range.par_mapv_inplace(|item| make_iteration(item, order));

    let p = range.map(|a| {
        if a.im == 0.0 || a.im == -0.0 {
            Complex64::new(a.re, 0.0)
        } else {
            *a
        }
    });
    let k = 1.0;

    Zpk { z, p, k }
}