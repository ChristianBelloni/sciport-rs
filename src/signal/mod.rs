use std::borrow::Cow;

use self::{
    band_filter::{lp2bf_zpk, BandFilter},
    output_type::{DesiredFilterOutput, FilterOutput, Zpk},
    tools::bilinear_zpk,
};

// filter design

pub mod band_filter;
/// Bessel/Thomson digital and analog filter design.
pub mod bessel;
/// Butterworth digital and analog filter design.
pub mod butter;
pub mod butterord;
pub mod cheby1;
pub mod cheby2;
pub mod ellip;
// end filter design
//
pub mod output_type;
pub mod tools;

#[derive(Debug, Clone, Copy)]
pub enum Analog {
    True,
    False { fs: f64 },
}

impl Analog {
    pub fn is_analog(&self) -> bool {
        match self {
            Self::True => true,
            Self::False { .. } => false,
        }
    }
}

/// Generic iir_filter
///
/// Takes a filter prototype and returns the final filter in the desired output
pub fn iir_filter(
    proto: Zpk,
    _order: u32,
    mut band_filter: BandFilter,
    mut analog: Analog,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let mut warped = band_filter;
    match &mut analog {
        Analog::True => {}
        Analog::False { fs } => {
            band_filter = (band_filter * 2.0) / *fs;
            *fs = 2.0;
            warped = ((band_filter * std::f64::consts::PI) / *fs).tan() * 2.0 * *fs;
        }
    }

    let mut result = lp2bf_zpk(proto, warped);
    if let Analog::False { fs } = &analog {
        result = bilinear_zpk(result, *fs);
    }

    FilterOutput::get_output(result, desired_output)
}

pub trait IIRFilter<T: ToOwned> {
    fn design_filter(&self) -> T;
    fn cache(&self) -> &Option<T>;
    fn cache_mut(&mut self) -> &mut Option<T>;

    fn get_filter(&mut self) -> &T {
        if self.cache().is_some() {
            return self.cache().as_ref().unwrap();
        } else {
            let filter = self.design_filter();
            self.cache_mut().insert(filter)
        }
    }
}

#[macro_export]
macro_rules! impl_iir {
    ($f:ty, $T:ty, $sel:ident, $design:expr) => {
        impl $crate::signal::IIRFilter<$T> for $f {
            fn cache(&self) -> &Option<$T> {
                &self.cache
            }

            fn cache_mut(&mut self) -> &mut Option<$T> {
                &mut self.cache
            }

            fn design_filter(&$sel) -> $T {
                $design
            }
        }
    };
}
