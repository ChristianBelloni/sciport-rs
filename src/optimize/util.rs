use num::complex::{Complex32, Complex64, ComplexFloat};

use crate::optimize::*;

pub fn approx_derivative<F, C, M>(fun: F) -> impl Fn(C) -> C
where
    F: Fn(C) -> C,
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
{
    let delta: C = C::epsilon().powf(C::from(1.0 / 3.0).unwrap().re());
    move |x: C| (fun(x + delta) - fun(x)) / delta
}
pub fn approx_second_derivative<F, C, M>(fun: F) -> impl Fn(C) -> C
where
    F: Fn(C) -> C,
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
{
    let delta: C = C::epsilon().powf(C::from(1.0 / 3.0).unwrap().re());
    move |x: C| (fun(x + delta) - fun(x) * C::from(2.0).unwrap() + fun(x - delta)) / delta / delta
}

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
