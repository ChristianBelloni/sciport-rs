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

mod convolution;
mod filter_design;
#[allow(unused)]
mod fir_filter_design;
mod signal_tools;

//pub use convolution::*;
pub use filter_design::*;

pub use fir_filter_design::{firwin, windows, Firwin1Filter, GenericFIRFilterSettings, WindowType};
use ndarray::{Array, Array1, Dimension, Ix1};
use num::Complex;

pub mod band_filter;
pub mod output_type;
pub mod tools;

pub use band_filter::{BandFilter, GenericBandFilter};
pub use filter_design::GenericIIRFilterSettings;
pub use filter_design::IIRFilterDesign;

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

pub trait Filter<T> {
    fn lfilter(
        &self,
        x: Array1<Complex<T>>,
        zi: Option<Array1<Complex<T>>>,
    ) -> LFilterOutput<T, Ix1>;
}

#[derive(Debug, Clone)]
pub struct LFilterOutput<T, D: Dimension> {
    pub filtered: Array<Complex<T>, D>,
    pub zi: Option<Array<Complex<T>, D>>,
}
