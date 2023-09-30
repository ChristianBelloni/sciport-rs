use num::complex::{Complex32, Complex64};

use crate::optimize::*;

/// A trait implemented for `f32`,`f64`,`Complex32` and `Complex64`
///
/// allow the type to be able to call `espilon()` for safe division
pub trait Espilon {
    fn epsilon() -> Self;
}
impl Espilon for f64 {
    fn epsilon() -> f64 {
        f64::EPSILON
    }
}
impl Espilon for f32 {
    fn epsilon() -> f32 {
        f32::EPSILON
    }
}
impl Espilon for Complex32 {
    fn epsilon() -> Complex32 {
        Complex32::new(f32::epsilon(), f32::epsilon())
    }
}
impl Espilon for Complex64 {
    fn epsilon() -> Complex64 {
        Complex64::new(f64::epsilon(), f64::epsilon())
    }
}
