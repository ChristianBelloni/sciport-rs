use crate::signal::band_filter::BandFilter;
use crate::signal::output_type::{Ba, DesiredFilterOutput, FilterOutput, Zpk};
use crate::signal::{iir_filter, Analog};
use ndarray::{array, Array1};
use num::complex::{Complex64, ComplexFloat};
use num::Complex;
use std::marker::PhantomData;

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
    iir_filter(dbg!(proto), order, band_filter, analog, desired_output)
}
pub fn cheb1ap(order: u32, rp: f64) -> Zpk {
    if order == 0 {
        let zpk = Zpk {
            z: Default::default(),
            p: Default::default(),
            k: dbg!(10.0_f64.powf(-rp / 20.0)),
        };
        return dbg!(zpk);
    }
    let rp = Complex::from(rp);
    let z: Array1<Complex64> = array![];
    let eps = Complex::from(10.0_f64.powc(0.1 * rp) - 1.0).sqrt();
    let mu = 1.0 / Complex::from(order as f64) * (Complex::from(1.0) / eps).asinh();
    let m: Array1<Complex64> = ((-(order as i32) + 1)..(order as i32))
        .step_by(2)
        .map(|a| (a as f64).into())
        .collect();

    let theta = m.map(|a| std::f64::consts::PI * a / (Complex::from(2.0) * order as f64));

    let p: Array1<Complex64> = theta.mapv(|a| -(Complex64::from(mu) + Complex64::i() * a).sinh());

    let mut k = p.map(|a| -a).product().re;

    if order % 2 == 0 {
        k /= (1.0 + eps.norm_sqr()).sqrt();
    }

    Zpk { z, p, k }
}

#[cfg(test)]
mod tests {

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