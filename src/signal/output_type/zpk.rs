use super::Ba;
use std::ops::Mul;

use num::{Complex, Num};

use super::{Filter, Zpk};

impl Zpk {}

impl<T: Num + Copy> Mul<T> for Zpk
where
    Complex<f64>: Mul<T>,
    T: Mul<Complex<f64>, Output = Complex<f64>>,
    T: Mul<f64, Output = f64>,
{
    type Output = Self;
    fn mul(mut self, rhs: T) -> Self::Output {
        self.z.iter_mut().for_each(|a| {
            *a = rhs * *a;
        });

        self.p.iter_mut().for_each(|a| {
            *a = rhs * *a;
        });

        self.k = rhs * self.k;

        self
    }
}

impl Filter for Zpk {
    fn lfilter<D: ndarray::Dimension>(
        &self,
        b: ndarray::Array1<f64>,
        a: ndarray::Array1<f64>,
        x: ndarray::Array<f64, D>,
        zi: Option<ndarray::Array<f64, D>>,
    ) -> Option<ndarray::Array<f64, D>> {
        let ba: Ba = self.clone().into();
        ba.lfilter(b, a, x, zi)
    }
}
