use ndarray::{array, Array1};
use num::{complex::ComplexFloat, traits::FloatConst, Complex, Float};

use crate::signal::output_type::GenericZpk;

use super::{GenericFilterSettings, ProtoFilter};

pub struct ButterFilter<T> {
    pub settings: GenericFilterSettings<T>,
}

impl<T: Float + FloatConst + ComplexFloat + Clone> ProtoFilter<T> for ButterFilter<T> {
    fn proto_filter(&self) -> crate::signal::output_type::GenericZpk<T> {
        buttap(self.settings.order)
    }

    fn filter_settings(&self) -> &GenericFilterSettings<T> {
        &self.settings
    }
}

pub fn buttap<T: Float>(order: u32) -> GenericZpk<T> {
    use std::f64::consts::PI;
    let order = order as i32;
    let z = array![];
    let range = Array1::range(
        T::from(1.0 - order as f64).unwrap(),
        T::from(order).unwrap(),
        T::from(2.0).unwrap(),
    );

    let range = range.mapv(|a| Complex::from(a));
    let k = T::one();

    let p = -(range
        .mapv(|a| Complex::i() * T::from(PI).unwrap() * a / T::from(2 * order).unwrap()))
    .mapv(Complex::exp);

    let p = p.mapv(|mut a| {
        if a.im == T::neg_zero() {
            a.im = T::zero()
        }
        a
    });

    GenericZpk { z, p, k }
}
