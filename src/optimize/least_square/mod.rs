use std::fmt::{Debug, Display};

use crate::odr::polynomial::{Polynomial, PolynomialCoef};
use crate::optimize::util::Espilon;
use nalgebra::{ComplexField, Scalar};
use num::complex::ComplexFloat;

/// # poly_fit
/// calculate the polynomial least square curve fit on data `x` and `y`
///
/// `order` define the order of the returned polynomial.
///
/// - if `order < x.len()`, the result will be least square curve fitting
///
/// - if `order = x.len()`, the result will be exact polynomial solve
///
/// - if `order > x.len()`, the result will be least sqaure of the coefficient of polynomial
///
/// ## Example
/// ```
/// # use sciport_rs::optimize::least_square::poly_fit;
///
/// let x = vec![1.0,2.0,3.0];
/// let y = vec![2.0,1.0,2.0];
///
/// let poly = poly_fit(&x,&y,2).unwrap();
/// ```
///
/// ## Errors
/// This function will return an error
/// - if the lenght if `x` and `y` are not the same.
/// - the svd solve fail.
pub fn poly_fit<'a, T, Q>(
    x: impl IntoIterator<Item = &'a T>,
    y: impl IntoIterator<Item = &'a T>,
    order: usize,
) -> Result<Polynomial<T>, String>
where
    T: Debug + Display + ComplexField<RealField = Q> + PolynomialCoef + Espilon,
    Q: ComplexFloat + Scalar + Debug + Espilon,
{
    let x = x.into_iter().collect::<Vec<&'a T>>();
    let y = y.into_iter().collect::<Vec<&'a T>>();

    if x.len() != y.len() {
        return Err(format!(
            "lsq_linear failed due to: len of x: {} and y: {} are not equal",
            x.len(),
            y.len()
        ));
    }

    let y = nalgebra::DVector::from_iterator(y.len(), y.into_iter().cloned());

    let rows = (0..x.len())
        .map(move |i| {
            nalgebra::DVector::from_iterator(
                order + 1,
                (0..(order + 1)).map(|a| ComplexFloat::powi(*x[i], a as i32)),
            )
            .transpose()
        })
        .collect::<Vec<_>>();

    nalgebra::DMatrix::from_rows(rows.as_slice())
        .svd(true, true)
        .solve(&y, Q::epsilon())
        .map(|s| s.into_iter().cloned().collect())
        .map_err(|e| format!("lsq_linear failed due to: {}", e))
}
