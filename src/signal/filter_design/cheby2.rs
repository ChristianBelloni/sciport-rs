use crate::signal::output_type::GenericZpk;
use ndarray::{array, concatenate, Array1, Axis};
use num::{complex::ComplexFloat, Complex, Float};

use super::{error::Infallible, GenericIIRFilterSettings, ProtoIIRFilter};

pub struct Cheby2Filter<T> {
    pub rs: T,
    pub settings: GenericIIRFilterSettings<T>,
}

impl<T: Float + num::traits::FloatConst + ComplexFloat + Clone> ProtoIIRFilter<T>
    for Cheby2Filter<T>
{
    fn proto_filter(&self) -> Result<GenericZpk<T>, crate::signal::error::Error> {
        Ok(cheb2ap(self.settings.order, self.rs).map_err(super::Error::from)?)
    }

    fn filter_settings(&self) -> &GenericIIRFilterSettings<T> {
        &self.settings
    }
}

pub fn cheb2ap<T: Float>(order: u32, rs: T) -> Result<GenericZpk<T>, Infallible> {
    use std::f64::consts::PI;
    if order == 0 {
        return Ok(GenericZpk {
            z: array![],
            p: array![],
            k: T::one(),
        });
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

    let k = match z.len().cmp(&p.len()) {
        std::cmp::Ordering::Less => {
            let tmp_z = concatenate![Axis(0), -&z, Array1::ones(p.len() - z.len())];
            (-&p / tmp_z).product()
        }
        std::cmp::Ordering::Equal => (-&p / -&z).product(),
        std::cmp::Ordering::Greater => {
            let tmp_p = concatenate![Axis(0), -&p, Array1::ones(z.len() - p.len())];
            (tmp_p / -&z).product()
        }
    }
    .re;
    Ok(GenericZpk { z, p, k })
}
