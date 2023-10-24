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
    band_filter::{lp2bf_zpk, BandFilter, GenericBandFilter},
    output_type::{DesiredFilterOutput, FilterGetOutputBounds, GenericFilterOutput, GenericZpk},
    tools::bilinear_zpk,
};

mod convolution;
mod filter_design;
mod fir_filter_design;
mod signal_tools;

//pub use convolution::*;
pub use filter_design::{bessel, butter, cheby1, cheby2};

pub use fir_filter_design::*;
use trait_set::trait_set;

pub mod band_filter;
pub mod output_type;
pub mod tools;

pub use filter_design::FilterDesign;
pub use filter_design::GenericFilterSettings;

pub type Sampling = GenericSampling<f64>;

#[derive(Debug, Clone, Copy)]
pub enum GenericSampling<T> {
    Analog,
    Digital { fs: T },
}

impl<T> GenericSampling<T> {
    pub fn is_analog(&self) -> bool {
        match self {
            Self::Analog => true,
            Self::Digital { .. } => false,
        }
    }

    pub fn cast<K>(self) -> GenericSampling<K>
    where
        K: From<T>,
    {
        match self {
            GenericSampling::Analog => GenericSampling::Analog,
            GenericSampling::Digital { fs } => GenericSampling::Digital { fs: fs.into() },
        }
    }

    pub fn cast_with_fn<K>(self, f: impl Fn(T) -> K) -> GenericSampling<K> {
        match self {
            GenericSampling::Analog => GenericSampling::Analog,
            GenericSampling::Digital { fs } => GenericSampling::Digital { fs: f(fs) },
        }
    }
}

trait_set! {
    pub trait IIRFilterBounds = FilterGetOutputBounds;
}

/// Generic iir_filter
///
/// Takes a filter prototype and returns the final filter in the desired output
pub fn iir_filter<T>(
    proto: GenericZpk<T>,
    _order: u32,
    mut band_filter: GenericBandFilter<T>,
    mut analog: GenericSampling<T>,
    desired_output: DesiredFilterOutput,
) -> GenericFilterOutput<T>
where
    T: IIRFilterBounds,
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

    GenericFilterOutput::get_output(result, desired_output)
}
