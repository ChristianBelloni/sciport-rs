use crate::{
    signal::{signal_tools::linear_filter, tools::zpk2ba},
    tools::convolve1d,
};
use ndarray::{Array1, Ix1};
use num::{complex::ComplexFloat, traits::FloatConst, Complex, Float, Zero};

use super::{Filter, GenericBa, GenericZpk, LFilterOutput};

impl<T: Float + FloatConst + ComplexFloat> From<GenericZpk<T>> for GenericBa<T> {
    fn from(value: GenericZpk<T>) -> Self {
        zpk2ba(value)
    }
}
#[allow(unused)]
fn mul_by_x(input: &mut Vec<Complex<f64>>) {
    input.push(Complex::zero());
}

#[allow(unused)]
fn mul_by_scalar(input: &mut [Complex<f64>], scalar: Complex<f64>) {
    input.iter_mut().for_each(|a| *a *= scalar);
}

#[allow(unused)]
fn sub(input: &mut [Complex<f64>], other: &[Complex<f64>]) {
    for (i, item) in other.iter().enumerate() {
        *input.get_mut(i).unwrap() = *input.get(i).unwrap() - item;
    }
}

impl<T: Float> Filter<T> for GenericBa<T> {
    fn lfilter(
        &self,
        x: ndarray::Array1<Complex<T>>,
        zi: Option<ndarray::Array1<Complex<T>>>,
    ) -> LFilterOutput<T, Ix1> {
        let b = self.b.clone();
        let a = self.a.clone();

        let zi = if let Some(zi) = zi {
            zi
        } else {
            Array1::zeros(b.raw_dim() - 1)
        };

        if a.len() == 1 {
            let b = b.mapv(|e| e / a[0]);

            let out_full = convolve1d(x.view(), b.view());

            LFilterOutput {
                filtered: out_full,
                zi: None,
            }
        } else {
            LFilterOutput {
                filtered: linear_filter(b, a, x, zi),
                zi: None,
            }
        }
    }
}
