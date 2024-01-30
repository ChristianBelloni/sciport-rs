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
        self.z = self.z.mapv(|a| rhs * a);
        self.p = self.p.mapv(|a| rhs * a);
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

#[cfg(test)]
mod tests {
    use ndarray::Array1;
    use ndarray_rand::{rand_distr::Normal, RandomExt};
    use num::Complex;
    use rand::{thread_rng, Rng};

    use crate::signal::output_type::GenericZpk;

    #[test]
    fn test_zpk_mul() {
        let z = Array1::random(6, Normal::new(10.0, 1.0).unwrap()).mapv(Complex::from);
        let p = Array1::random(6, Normal::new(10.0, 1.0).unwrap()).mapv(Complex::from);
        let k: f64 = thread_rng().gen_range(0.1..1.2);

        let mut zpk = GenericZpk { z, p, k };
        let rhs = 6.0;
        let mul_zpk = zpk.clone() * rhs;

        zpk.z = zpk.z.mapv(|a| rhs * a);
        zpk.p = zpk.p.mapv(|a| rhs * a);
        zpk.k *= rhs;

        assert_eq!(mul_zpk, zpk);
    }
}
