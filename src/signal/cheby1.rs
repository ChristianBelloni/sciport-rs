use std::marker::PhantomData;

use super::band_filter::BandFilter;
use super::output_type::{Ba, DesiredFilterOutput, FilterOutput, Zpk};
use super::{iir_filter, Analog};
use num::complex::Complex64;
use num::One;

pub struct Cheby1FilterStandalone<T>(PhantomData<T>);

impl Cheby1FilterStandalone<Zpk> {
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog, rp: f64) -> Zpk {
        cheby1_filter(order, band_filter, analog, rp, DesiredFilterOutput::Zpk).zpk()
    }
}

impl Cheby1FilterStandalone<Ba> {
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog, rp: f64) -> Ba {
        cheby1_filter(order, band_filter, analog, rp, DesiredFilterOutput::Ba).ba()
    }
}

pub struct Cheby1Filter<T> {
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    rp: f64,
    cache: Option<T>,
}

crate::impl_iir!(
    Cheby1Filter::<Zpk>,
    Zpk,
    self,
    Cheby1FilterStandalone::<Zpk>::filter(self.order, self.band_filter, self.analog, self.rp)
);

crate::impl_iir!(
    Cheby1Filter::<Ba>,
    Ba,
    self,
    Cheby1FilterStandalone::<Ba>::filter(self.order, self.band_filter, self.analog, self.rp)
);

pub(crate) fn cheby1_filter(
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    rp: f64,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let proto = cheb1ap(order, rp);
    iir_filter(proto, order, band_filter, analog, desired_output)
}
pub fn cheb1ap(order: u32, rp: f64) -> Zpk {
    if order == 0 {
        let zpk = Zpk {
            z: Default::default(),
            p: Default::default(),
            k: 10.0_f64.powf(-rp / 20.0),
        };
        return zpk;
    }

    let z: Vec<Complex64> = vec![];
    let eps = (10.0_f64.powf(0.1 * rp) - 1.0).sqrt();
    let mu = 1.0 / (order as f64) * (1.0 / eps).asinh();
    let m = ((-(order as i32) + 1)..(order as i32)).step_by(2);

    let theta = m.map(|a| std::f64::consts::PI * (a as f64) / (2.0 * order as f64));

    let p: Vec<_> = theta.map(|a| -(mu + Complex64::i() * a).sinh()).collect();

    let mut k = p
        .iter()
        .map(|&a| -a)
        .fold(Complex64::one(), |acc, item| acc * item)
        .re;

    if order % 2 == 0 {
        k /= (1.0 + eps * eps).sqrt();
    }

    Zpk { z, p, k }
}

#[cfg(test)]
mod tests {
    use num::Signed;
    macro_rules! assert_almost_vec_eq {
        ($v1: expr, $v2: expr, $tol: expr) => {
            for (i, item) in $v1.iter().enumerate() {
                let err = (item - $v2[i]).norm().abs();
                println!("{err}");
                let is_in_tol = err < $tol;
                assert!(is_in_tol);
            }
        };
    }
    use super::*;
    #[test]
    fn test_cheb1ap() {
        // [-0.17221491+0.9159552j, -0.34442981-0.j,-0.17221491-0.9159552j]
        let res = cheb1ap(3, 2.3);

        assert_almost_vec_eq!(
            &res.p,
            &[
                Complex64::new(-0.17221491, 0.9159552),
                Complex64::new(-0.34442981, -0.0),
                Complex64::new(-0.17221491, -0.9159552)
            ],
            10.0_f64.powi(-8)
        )
    }
}
