use super::{GenericBa, GenericZpk, LFilterOutput};
use std::ops::Mul;

use ndarray::Ix1;
use num::{traits::FloatConst, Complex, Float, Num};

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

impl<T: Float + FloatConst> Filter<T> for GenericZpk<T> {
    fn lfilter(
        &self,
        x: ndarray::Array1<Complex<T>>,
        zi: Option<ndarray::Array1<Complex<T>>>,
    ) -> LFilterOutput<T, Ix1> {
        let ba: GenericBa<T> = self.clone().into();
        ba.lfilter(x, zi)
    }
}
