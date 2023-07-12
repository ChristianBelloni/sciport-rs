use std::{borrow::Cow, marker::PhantomData};

use num::Complex;

use super::{
    band_filter::BandFilter,
    iir_filter,
    output_type::{Ba, DesiredFilterOutput, FilterOutput, Zpk},
    Analog, IIRFilter,
};

pub(crate) fn butter_filter(
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let proto = buttap(order);
    iir_filter(proto, order, band_filter, analog, desired_output)
}
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

pub fn buttap(order: u32) -> Zpk {
    let order = order as i32;
    let z = vec![];
    let range = ((-order + 1)..order).step_by(2);

    fn make_iteration(item: i32, order: i32) -> Complex<f64> {
        let numerator = Complex::new(0.0, 1.0) * std::f64::consts::PI * (item as f64);
        let denominator = 2.0 * order as f64;

        let temp = numerator / denominator;

        -temp.exp()
    }

    let p = range
        .map(|item| make_iteration(item, order))
        .collect::<Vec<_>>();

    let k = 1.0;

    Zpk { z, p, k }
}

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
