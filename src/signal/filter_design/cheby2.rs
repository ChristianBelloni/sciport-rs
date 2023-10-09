use crate::signal::band_filter::BandFilter;
use crate::signal::output_type::{Ba, DesiredFilterOutput, FilterOutput, Zpk};
use crate::signal::{iir_filter, Analog};
use ndarray::{array, Array1};
use num::complex::Complex64;
use std::marker::PhantomData;

pub struct Cheby2FilterStandalone<T>(PhantomData<T>);

impl Cheby2FilterStandalone<Zpk> {
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog, rs: f64) -> Zpk {
        cheby2_filter(order, band_filter, analog, rs, DesiredFilterOutput::Zpk).zpk()
    }
}

impl Cheby2FilterStandalone<Ba> {
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog, rs: f64) -> Ba {
        cheby2_filter(order, band_filter, analog, rs, DesiredFilterOutput::Ba).ba()
    }
}

pub struct Cheby2Filter<T> {
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    rs: f64,
    cache: Option<T>,
}

crate::impl_iir!(
    Cheby2Filter::<Zpk>,
    Zpk,
    self,
    Cheby2FilterStandalone::<Zpk>::filter(self.order, self.band_filter, self.analog, self.rs)
);

crate::impl_iir!(
    Cheby2Filter::<Ba>,
    Ba,
    self,
    Cheby2FilterStandalone::<Ba>::filter(self.order, self.band_filter, self.analog, self.rs)
);

pub(crate) fn cheby2_filter(
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    rs: f64,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let proto = cheb2ap(order, rs);
    iir_filter(proto, order, band_filter, analog, desired_output)
}

pub fn cheb2ap(order: u32, rs: f64) -> Zpk {
    if order == 0 {
        return Zpk {
            z: array![],
            p: array![],
            k: 1.0,
        };
    }

    let de = 2.0 / (10.0_f64.powf(0.1 * rs) - 1.0).sqrt();

    let mu = (2.0 / de).asinh() / order as f64;
    let m: Array1<Complex64> = if order % 2 == 0 {
        let first = ((-(order as i32) + 2)..0)
            .step_by(2)
            .map(|a| (a as f64).into());
        let second = (2..(order as i32)).step_by(2).map(|a| (a as f64).into());

        first.chain(second).collect()
    } else {
        ((-(order as i32) + 2)..(order as i32))
            .step_by(2)
            .map(|a| (a as f64).into())
            .collect()
    };

    let z = m.map(|a| {
        -(Complex64::i() / (a * std::f64::consts::PI / (2.0 * (order as f64))).sin()).conj()
    });

    let tmp: Array1<Complex64> = ((-(order as i32) + 2)..(order as i32))
        .map(|a| (a as f64).into())
        .collect();

    let mut p =
        tmp.map(|a| -(Complex64::i() * std::f64::consts::PI * a / (2.0 * order as f64)).exp());

    p.map_inplace(|a| {
        *a = mu.sinh() * a.re + Complex64::i() * mu.cosh() * a.im;
    });

    let k = (p.map(|a| -a).product() / z.map(|a| -a).product()).re;

    Zpk { z, p, k }
}
