//! # Signal Processing
//! The signal processing toolbox currently contains some filtering functions, a limited set of filter design tools,<br/>
//! and a few B-spline interpolation algorithms for 1- and 2-D data. While the B-spline algorithms could<br/>
//! technically be placed under the interpolation category, they are included here because they only work with<br/>
//! equally-spaced data and make heavy use of filter-theory and transfer-function formalism to provide a fast B-spline transform.<br/>
//! To understand this section, you will need to understand that a signal in sciport-rs is an array of real or complex numbers.
//!
//! ## Filter Design
//!
//! Time-discrete filters can be classified into finite response (FIR) filters and infinite response (IIR) filters.<br/>
//! FIR filters can provide a linear phase response, whereas IIR filters cannot.<br/>
//! sciport-rs provides functions for designing both types of filters.
//!
//! ### IIR Filter
//!
//! sciport-rs provides two functions to directly design IIR iirdesign and iirfilter, where the filter type (e.g., elliptic)<br/>
//! is passed as an argument and several more filter design functions for specific filter types, e.g., ellip.
//! ### Filter coefficients
//!
//! Filter coefficients can be stored in several different formats:
//!  - [`Ba`](`crate::signal::output_type::Ba`)
//!  - [`Zpk`](`crate::signal::output_type::Zpk`)
//!  - [`Sos`](`crate::signal::output_type::Sos`) (currently unsupported)
//!
//! # References:
//!
//! The documentation on this page is largely been copied from the [SciPy](https://docs.scipy.org/doc/scipy/tutorial/signal.html) documentation
use self::{
    band_filter::{lp2bf_zpk, BandFilter},
    output_type::{DesiredFilterOutput, FilterOutput, Zpk},
    tools::bilinear_zpk,
};

mod convolution;
mod filter_design;

//pub use convolution::*;
pub use filter_design::*;

pub mod band_filter;
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

    /// Computes and returns the filter coefficients
    fn get_filter(&mut self) -> &T {
        if self.cache().is_some() {
            return self.cache().as_ref().unwrap();
        } else {
            let filter = self.design_filter();
            self.cache_mut().insert(filter)
        }
    }
}
#[doc(hidden)]
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
