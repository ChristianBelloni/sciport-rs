use crate::signal::output_type::GenericZpk;
use crate::tools::complex::normalize_zeros;
use ndarray::{array, Array1};
use num::complex::ComplexFloat;
use num::*;
use std::f64::consts::PI;

use super::{GenericFilterSettings, ProtoFilter};

pub struct Cheby1Filter<T> {
    pub rp: T,
    pub settings: GenericFilterSettings<T>,
}

impl<T: Float + traits::FloatConst + ComplexFloat + Clone> ProtoFilter<T> for Cheby1Filter<T> {
    fn proto_filter(&self) -> crate::signal::output_type::GenericZpk<T> {
        cheb1ap(self.settings.order, self.rp)
    }

    fn filter_settings(&self) -> &GenericFilterSettings<T> {
        &self.settings
    }
}

pub fn cheb1ap<T>(order: u32, rp: T) -> GenericZpk<T>
where
    T: Float,
{
    let from = |v: f64| -> T { NumCast::from(v).unwrap() };

    if order == 0 {
        return GenericZpk {
            z: array![],
            p: array![],
            k: from(10.0).powf(-rp / from(20.0)),
        };
    }

    let z = array![];

    let eps = (from(10.0).powf(from(0.1) * rp) - from(1.0)).sqrt();
    let mu = from(1.0) / from(order as _) * (from(1.0) / eps).asinh();

    let m = Array1::range(from(1.0 - (order as f64)), from(order as _), from(2.0));
    let theta = m
        .mapv(|a| a * from(PI))
        .mapv(|a| a / from(2.0 * (order as f64)));
    let p = -(theta.mapv(|a| Complex::new(mu, from(0.0)) + Complex::i() * a)).mapv(Complex::sinh);

    let mut k = (-&p).product().re;
    if order % 2 == 0 {
        k = k / (from(1.0) + eps * eps).sqrt();
    }
    let z = normalize_zeros(z);
    let p = normalize_zeros(p);

    GenericZpk { z, p, k }
}
