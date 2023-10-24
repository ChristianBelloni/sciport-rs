use num::{complex::ComplexFloat, traits::FloatConst, Float};

use super::{
    band_filter::GenericBandFilter,
    iir_filter,
    output_type::{DesiredFilterOutput, GenericFilterOutput, GenericZpk},
    GenericSampling,
};

pub mod bessel;
pub mod butter;
// pub mod butterord;
pub mod cheby1;
pub mod cheby2;
// pub mod ellip;

pub use bessel::*;
pub use butter::*;
pub use cheby1::*;
pub use cheby2::*;

pub struct GenericIIRFilterSettings<T> {
    pub order: u32,
    pub band_filter: GenericBandFilter<T>,
    pub analog: GenericSampling<T>,
}

pub trait ProtoIIRFilter<T: Float + FloatConst + ComplexFloat> {
    fn proto_filter(&self) -> GenericZpk<T>;

    fn filter_settings(&self) -> &GenericIIRFilterSettings<T>;
}

pub trait IIRFilterDesign<T: Float + FloatConst + ComplexFloat>: ProtoIIRFilter<T> {
    fn compute_filter(&self, desired_output: DesiredFilterOutput) -> GenericFilterOutput<T> {
        let proto = self.proto_filter();
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
