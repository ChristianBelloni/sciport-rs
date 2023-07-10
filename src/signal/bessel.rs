use std::marker::PhantomData;

use num::{complex::ComplexFloat, Complex, One, Zero};

use crate::{signal::tools::polyval, special::kv::kve};

use super::{
    band_filter::BandFilter,
    iir_filter,
    output_type::{DesiredFilterOutput, FilterOutput, Zpk},
    Analog,
};

#[derive(Debug, Clone, Copy)]
pub enum BesselNorm {
    Phase,
    Delay,
    Mag,
}

pub(crate) fn bessel_filter(
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    norm: BesselNorm,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let proto = besselap(order, norm);

    iir_filter(proto, order, band_filter, analog, desired_output)
}

pub struct BesselFilter<T>(PhantomData<T>);

impl BesselFilter<Zpk> {
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog, norm: BesselNorm) -> Zpk {
        bessel_filter(order, band_filter, analog, norm, DesiredFilterOutput::Zpk).zpk()
    }
}

fn besselap(order: u32, norm: BesselNorm) -> Zpk {
    let z = Vec::<Complex<f64>>::new();
    let mut p: Vec<Complex<f64>>;
    let k = 1.0;
    if order == 0 {
        p = vec![];
    } else {
        println!("making zeroes");
        p = _bessel_zeros(order).into_iter().map(|a| 1.0 / a).collect();
        println!("making factorials");
        let a_last = _falling_factorial(2 * order, order);

        p.iter_mut()
            .for_each(|a| *a = Complex::from(10.0).powf(-(a_last.log10() / order as f64)));
    }

    Zpk { z, p, k }
}

fn _falling_factorial(x: u32, n: u32) -> f64 {
    let mut y = 1.0;

    for i in (x - n + 1)..(x + 1) {
        y *= i as f64;
    }
    y
}

fn _bessel_zeros(order: u32) -> Vec<Complex<f64>> {
    if order == 0 {
        return Default::default();
    }

    let x0 = _campos_zeros(order);
    println!("made campos_zeroes");
    let f = |x: Complex<f64>| kve(order as f64 - 0.5, x);

    let fp = |x: Complex<f64>| {
        let order = order as f64;
        (kve(order - 0.5, 1.0 / x) / (2.0 * x.powi(2))) - (kve(order + 0.5, 1.0 / x) / (x.powi(2)))
            + (kve(order + 1.5, 1.0 / x) / (2.0 * x.powi(2)))
    };
    let x = _aberth(f, fp, &x0);
    // TODO implement newton's method to improve accuracy and average complex conjugates

    x
}

fn _aberth<F: Fn(Complex<f64>) -> Complex<f64>, FP: Fn(Complex<f64>) -> Complex<f64>>(
    f: F,
    fp: FP,
    x0: &[Complex<f64>],
) -> Vec<Complex<f64>> {
    let tol = 10.0_f64.powi(-15);
    let maximiter = 50;

    let n = x0.len();
    let mut beta = vec![Complex::<f64>::zero(); n];

    let mut x = x0.to_vec();
    println!("{x:#?}");

    for iteration in 0..maximiter {
        println!("iteration {iteration} of {maximiter}");
        let alpha = x0.iter().map(|&a| -f(a) / fp(a)).collect::<Vec<_>>();
        for i in 0..n {
            beta[i] = Complex::from(1.0_f64)
                / x[(i + 1)..]
                    .iter()
                    .map(|&a| -a + x[i])
                    .sum::<Complex<f64>>();
            beta[i] +=
                Complex::from(1.0_f64) / x[..i].iter().map(|&a| -a + x[i]).sum::<Complex<f64>>();

            x.iter_mut()
                .enumerate()
                .for_each(|(i, a)| *a += alpha[i] / (1.0 + alpha[i] + beta[i]));

            if alpha.iter().all(|&a| a.abs() <= tol) {
                break;
            }
        }
    }
    println!("{x:#?}");
    x
}

fn _campos_zeros(order: u32) -> Vec<Complex<f64>> {
    if order == 0 {
        return vec![-Complex::one()];
    }
    let n = order as _;

    let s = polyval(n, [0.0, 0.0, 2.0, 0.0, -3.0, 1.0]);
    let b3 = polyval(n, [16.0, -8.0]) / s;
    let b2 = polyval(n, [-24.0, -12.0, 12.0]) / s;
    let b1 = polyval(n, [8.0, 24.0, -12.0, -2.0]) / s;
    let b0 = polyval(n, [0.0, -6.0, 0.0, 5.0, -1.0]) / s;

    let r = polyval(n, [0.0, 0.0, 2.0, 1.0]);

    let a1 = polyval(n, [-6.0, -6.0]) / r;
    let a2 = 6.0 / r;

    let k = 1..(order + 1);
    let x = k
        .clone()
        .map(|a| polyval(Complex::from(a as f64), [0.0.into(), a1, a2]))
        .collect::<Vec<_>>();
    let y = k
        .map(|a| polyval(Complex::from(a as f64), [b0, b1, b2, b3]))
        .collect::<Vec<_>>();
    println!("x: {x:#?}");
    println!("y: {y:#?}");
    assert_eq!(x.len(), y.len());
    x.iter()
        .zip(y)
        .map(|(x, y)| *x + Complex::new(0.0, 1.0) * y)
        .collect::<Vec<_>>()
}

#[cfg(test)]
mod tests {
    use crate::signal::band_filter::BandFilter;

    use super::{besselap, BesselFilter, BesselNorm};

    #[test]
    fn test_bessel_filter() {
        let filter = BesselFilter::filter(
            4,
            BandFilter::Lowpass(0.2),
            crate::signal::Analog::False { fs: 2.0 },
            BesselNorm::Phase,
        );

        

        println!("{filter:#?}");
    }

    #[test]
    fn test_besselap() {
        let res = besselap(4, BesselNorm::Phase);

        println!("{res:#?}");
    }
}
