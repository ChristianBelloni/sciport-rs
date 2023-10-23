use crate::odr::polynomial::{Polynomial, PolynomialCoef};
use crate::optimize::root_scalar::*;
use crate::optimize::util::Espilon;
use num::complex::{Complex32, Complex64, ComplexFloat};

pub trait IntoComplex<C>: ComplexFloat {
    fn into_complex(&self) -> C;
}
impl IntoComplex<Complex32> for f32 {
    fn into_complex(&self) -> Complex32 {
        Complex32::new(*self, 0.0)
    }
}
impl IntoComplex<Complex64> for f64 {
    fn into_complex(&self) -> Complex64 {
        Complex64::new(*self, 0.0)
    }
}
impl IntoComplex<Complex32> for Complex32 {
    fn into_complex(&self) -> Complex32 {
        *self
    }
}
impl IntoComplex<Complex64> for Complex64 {
    fn into_complex(&self) -> Complex64 {
        *self
    }
}
/// Solve quadratic equation for given equation
///
/// for `f32`, `f64`, it will return None if the determiniate is non-square-root-able
///
/// for `Complex32`, `Complex64`, will always return root.
pub fn quadratic_root<T>(a: T, b: T, c: T) -> Option<Vec<T>>
where
    T: ComplexFloat,
{
    let d = (b.powi(2) - a * c * T::from(4.0).unwrap()).sqrt();
    if d.is_nan() {
        return None;
    }
    let a2 = a * T::from(2.0).unwrap();
    Some(vec![(-b + d) / a2, (-b - d) / a2])
}

/// solve for all root of a given polynomial.
///
/// it always convert the given polynomial to complex domain first,
/// thus there will always be root
///
/// this function use halley's method to solve for roots
/// and deflate polynomial recursively,
/// it might panic if halley's method cannot solve for the root within 1000 iterations.
pub fn polynomial_roots<T, C, M>(polynomial: &Polynomial<T>) -> Vec<C>
where
    T: PolynomialCoef + Espilon + IntoMetric<M> + IntoComplex<C>,
    C: PolynomialCoef + Espilon + IntoMetric<M> + IntoComplex<C>,
    M: Metric,
{
    let polynomial: Polynomial<C> = polynomial
        .iter()
        .map(|&c| c.into_complex())
        .collect::<Polynomial<C>>()
        .saturate();

    match polynomial.degree() {
        0 => vec![],
        1 => vec![-polynomial[0] / polynomial[1]],
        2 => quadratic_root(polynomial[2], polynomial[1], polynomial[0]).unwrap(),
        _ => {
            let criteria = Some(
                OptimizeCriteria::empty()
                    .set_fltol(Some(M::from(1e-15).unwrap()))
                    .set_maxiter(Some(1000)),
            );

            let x0 = C::zero();

            let fun = polynomial.clone().as_fn();
            let dfun = polynomial.clone().differentiate().as_fn();
            let ddfun = polynomial.clone().differentiate().differentiate().as_fn();

            let res = halley_method(fun, dfun, ddfun, x0, criteria);

            let root = res.sol_x.unwrap();

            let deflated = polynomial.deflate(root).unwrap().0.to_owned();

            let roots: Vec<C> = polynomial_roots(&deflated);

            std::iter::once(root).chain(roots).collect()
        }
    }
}
