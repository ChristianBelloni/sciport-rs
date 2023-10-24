pub mod firwin1;
mod pass_zero;
mod tools;
pub mod windows;

pub use tools::*;

pub use self::windows::WindowType;
use super::{band_filter::GenericBandFilter, GenericSampling};

pub struct GenericFIRFilterSettings<T> {
    pub numtaps: i64,
    pub cutoff: GenericBandFilter<T>,
    pub width: Option<f64>,
    pub window: WindowType<T>,
    pub scale: bool,
    pub sampling: GenericSampling<T>,
}

pub use firwin1::*;
