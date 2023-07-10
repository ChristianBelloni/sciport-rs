use std::ops::Mul;

use num::{Complex, Num};

use super::Zpk;

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
