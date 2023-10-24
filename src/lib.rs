//! # Sciport-rs
//!
//! Sciport is a collection of mathematical algorithms ported from the popular python package Scipy
//!
//! # Api design
//!
//! The main philosophy behind sciport is to change the api surface of scipy to better utilize the
//! rich rust typesystem, when deciding between keeping the original function signature and
//! rewriting it to better represent the valid input space, more often than not we'll decide to
//! change it.<br/>
//! for example this is the scipy butter filter api:
//!
//! ```python
//! scipy.signal.butter(N: int, Wn: array_like, btype: String, analog: bool, output: String, fs:
//! float)
//! ```
//!
//! Wn represents a single or a pair of frequencies and btype is the type of filter,
//! however, a single frequency makes sense only for a subset of btypes and so does a pair,
//! in our implementation we rewrite this function like:
//!
//! ```
//! # use sciport_rs::signal::Sampling;
//! # use sciport_rs::signal::band_filter::BandFilter;
//! fn filter<T>(order: u32, band_filter: BandFilter, analog: Sampling) {  }
//! ```
//!
//! where T represents the output representation of the filter (Zpk, Ba, Sos), band_filter
//! encapsulates the original Wn and btype like this:
//!
//! ```
//!
//! pub enum BandFilter {
//!     Highpass(f64),
//!     Lowpass(f64),
//!     Bandpass { low: f64, high: f64 },
//!     Bandstop { low: f64, high: f64 },
//! }
//! ```
//!
//! and analog encapsulates analog and fs (since a sampling rate makes sense only when talking
//! about a digital filter) like this:
//!
//! ```
//! pub enum Sampling {
//!     Analog,
//!     Digital {
//!         fs: f64
//!     }
//! }
//! ```
//!
//! # Modules
//!
//! ## Signal Processing
//!
//! The signal processing toolbox currently contains some filtering functions, a limited set of filter design tools, and a few B-spline interpolation algorithms for 1- and 2-D data. While the B-spline algorithms could technically be placed under the interpolation category, they are included here because they only work with equally-spaced data and make heavy use of filter-theory and transfer-function formalism to provide a fast B-spline transform.
//!
//! ## Special
//!
//! The main feature of this module is the definition of numerous special functions
//! of mathematical physics. Available functions include airy, elliptic, bessel, gamma, beta,
//! hypergeometric, parabolic cylinder, mathieu, spheroidal wave, struve, and kelvin.
//!
#[allow(unused)]
pub mod odr;
#[allow(unused)]
pub mod optimize;
pub mod signal;
pub mod special;

pub(crate) mod tools;
