use super::{
    band_filter::GenericBandFilter,
    output_type::{DesiredFilterOutput, GenericFilterOutput, GenericZpk},
    GenericSampling,
};
use crate::signal::{band_filter::lp2bf_zpk, tools::bilinear_zpk};
use num::{complex::ComplexFloat, traits::FloatConst, Float};

pub mod bessel;
pub mod butter;
#[allow(unused)]
pub mod butterord;
pub mod cheby1;
pub mod cheby2;
pub mod error;
// pub mod ellip;
pub use error::Error;

pub use bessel::{besselap, BesselFilter, BesselNorm};
pub use butter::{buttap, ButterFilter};
pub use cheby1::{cheb1ap, Cheby1Filter};
pub use cheby2::{cheb2ap, Cheby2Filter};

/// Generic iir_filter
///
/// Takes a filter prototype and returns the final filter in the desired output
pub fn iir_filter<T>(
    proto: GenericZpk<T>,
    _order: u32,
    mut band_filter: GenericBandFilter<T>,
    mut analog: GenericSampling<T>,
    desired_output: DesiredFilterOutput,
) -> Result<GenericFilterOutput<T>, super::error::Error>
where
    T: Float + FloatConst + ComplexFloat,
{
    use std::f64::consts::PI;

    let mut warped: GenericBandFilter<T> = band_filter;
    match &mut analog {
        GenericSampling::Analog => {}
        GenericSampling::Digital { fs } => {
            band_filter = (band_filter * T::from(2).unwrap()) / *fs;
            *fs = T::from(2).unwrap();
            let tmp: GenericBandFilter<_> = ((band_filter * T::from(PI).unwrap()) / *fs).tan();
            warped = tmp * T::from(4).unwrap();
        }
    }

    let mut result = lp2bf_zpk(proto, warped);

    if let GenericSampling::Digital { fs } = &analog {
        result = bilinear_zpk(result, *fs);
    }

    Ok(GenericFilterOutput::get_output(result, desired_output))
}

pub struct GenericIIRFilterSettings<T> {
    pub order: u32,
    pub band_filter: GenericBandFilter<T>,
    pub analog: GenericSampling<T>,
}

pub trait ProtoIIRFilter<T: Float + FloatConst + ComplexFloat> {
    fn proto_filter(&self) -> Result<GenericZpk<T>, crate::signal::error::Error>;

    fn filter_settings(&self) -> &GenericIIRFilterSettings<T>;
}

pub trait IIRFilterDesign<T: Float + FloatConst + ComplexFloat>: ProtoIIRFilter<T> {
    fn compute_filter(
        &self,
        desired_output: DesiredFilterOutput,
    ) -> Result<GenericFilterOutput<T>, super::error::Error> {
        let proto = self.proto_filter()?;
        let settings = self.filter_settings();
        iir_filter(
            proto,
            settings.order,
            settings.band_filter,
            settings.analog,
            desired_output,
        )
    }
}

impl<T: Float + FloatConst + ComplexFloat, K> IIRFilterDesign<T> for K where K: ProtoIIRFilter<T> {}

pub trait OrdCompute<T> {
    fn compute_order(&self) -> Result<OrdResult<T>, crate::signal::error::Error>;
}

pub struct OrdResult<T> {
    pub order: u32,
    pub filter: GenericBandFilter<T>,
}
