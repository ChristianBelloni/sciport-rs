use crate::signal::output_type::GenericZpk;
use ndarray::{array, concatenate, Array1, Axis};
use num::{complex::ComplexFloat, Complex, Float};

use super::{GenericFilterSettings, ProtoFilter};

pub struct Cheby2Filter<T> {
    pub rs: T,
    pub settings: GenericFilterSettings<T>,
}

impl<T: Float + num::traits::FloatConst + ComplexFloat + Clone> ProtoFilter<T> for Cheby2Filter<T> {
    fn proto_filter(&self) -> crate::signal::output_type::GenericZpk<T> {
        cheb2ap(self.settings.order, self.rs)
    }

    fn filter_settings(&self) -> &GenericFilterSettings<T> {
        &self.settings
    }
}

pub fn cheb2ap<T: Float>(order: u32, rs: T) -> GenericZpk<T> {
    use std::f64::consts::PI;
    if order == 0 {
        return GenericZpk {
            z: array![],
            p: array![],
            k: T::one(),
        };
    }

    let de = T::one() / (T::from(10).unwrap().powf(T::from(0.1).unwrap() * rs) - T::one()).sqrt();
    let mu = (T::one() / de).asinh() / T::from(order).unwrap();

    let m = if order % 2 == 1 {
        concatenate![
            Axis(0),
            Array1::range(
                -T::from(order).unwrap() + T::one(),
                T::zero(),
                T::from(2).unwrap()
            ),
            Array1::range(
                T::from(2).unwrap(),
                T::from(order).unwrap(),
                T::from(2).unwrap()
            )
        ]
    } else {
        Array1::range(
            -T::from(order).unwrap() + T::one(),
            T::from(order).unwrap(),
            T::from(2).unwrap(),
        )
    }
    .mapv(Complex::from);
    let z = -m
        .mapv(|a| a * T::from(PI).unwrap())
        .mapv(|a| a / (T::from(2 * order).unwrap()))
        .mapv(Complex::sin)
        .mapv(|a| Complex::<T>::i() / a)
        .map(Complex::conj);

    let order = T::from(order).unwrap();

    let p = -Array1::range(-order + T::one(), order, T::from(2).unwrap())
        .mapv(|a| a / (T::from(2).unwrap() * order))
        .mapv(|a| Complex::i() * T::from(PI).unwrap() * a)
        .mapv(Complex::exp);
    let p = p.mapv(|a| {
        let re = mu.sinh() * a.re;
        let im = mu.cosh() * a.im;
        Complex::new(re, im)
    });
    let p = p.map(Complex::inv);

    let k = if z.len() == p.len() {
        (-&p / -&z).product()
    } else if z.len() > p.len() {
        let tmp_p = concatenate![Axis(0), -&p, Array1::ones(z.len() - p.len())];
        (tmp_p / -&z).product()
    } else {
        let tmp_z = concatenate![Axis(0), -&z, Array1::ones(p.len() - z.len())];
        (-&p / tmp_z).product()
    }
    .re;

    GenericZpk { z, p, k }
}
