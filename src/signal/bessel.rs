use std::marker::PhantomData;

use ndarray::{array, Array1};
use num::{
    complex::{Complex64, ComplexFloat},
    Complex, One, Zero,
};

use crate::{
    signal::tools::{newton, polyval},
    special::kve,
};

use super::{
    band_filter::BandFilter,
    iir_filter,
    output_type::{Ba, DesiredFilterOutput, FilterOutput, Zpk},
    Analog,
};

pub struct BesselFilter<T> {
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    norm: BesselNorm,
    cache: Option<T>,
}

crate::impl_iir!(
    BesselFilter<Zpk>,
    Zpk,
    self,
    BesselFilterStandalone::<Zpk>::filter(self.order, self.band_filter, self.analog, self.norm)
);

crate::impl_iir!(
    BesselFilter<Ba>,
    Ba,
    self,
    BesselFilterStandalone::<Ba>::filter(self.order, self.band_filter, self.analog, self.norm)
);

#[derive(Debug, Clone, Copy)]
pub enum BesselNorm {
    Phase,
    Delay,
    Mag,
}

/// Bessel/Thomson digital and analog filter design.
/// Design a Nth-order digital or analog Bessel filter and return the filter coefficients
///
/// # Notes
///
/// Also known as a Thomson filter, the analog Bessel filter has maximally flat group delay and<br/>
/// maximally linear phase response, with very little ringing in the step response.<br/>

/// The Bessel is inherently an analog filter. This function generates digital Bessel filters using the bilinear transform,<br/>
/// which does not preserve the phase response of the analog filter. As such, it is only approximately correct at frequencies below about fs/4.<br/>
/// To get maximally-flat group delay at higher frequencies, the analog Bessel filter must be transformed using phase-preserving techniques.<br/>
///
/// See [`besselap`] for implementation details and references.
///
/// # References
///
/// - Thomson, W.E., “Delay Networks having Maximally Flat Frequency Characteristics”, Proceedings</br>
/// of the Institution of Electrical Engineers, Part III, November 1949, Vol. 96, No. 44, pp. 487-490.
pub struct BesselFilterStandalone<T>(PhantomData<T>);

impl BesselFilterStandalone<Zpk> {
    /// Bessel/Thomson digital and analog filter design.
    /// Design a Nth-order digital or analog Bessel filter and return the filter coefficients
    ///
    /// # Parameters
    ///  - order : Order of the filter
    ///  - band: [`BandFilter`] to apply, for analog filters band is expressed as an angular</br>
    /// frequency (rad/s).
    /// for digital filters band is in the same units as [`Analog`]
    ///  - analog: Analog or Digital filter selection, when digital is selected<br/>
    ///  a sampling rate is required
    ///  - norm: Critical frequency normalization, see [`BesselNorm`] for more info
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog, norm: BesselNorm) -> Zpk {
        bessel_filter(order, band_filter, analog, norm, DesiredFilterOutput::Zpk).zpk()
    }
}

impl BesselFilterStandalone<Ba> {
    /// Bessel/Thomson digital and analog filter design.
    /// Design a Nth-order digital or analog Bessel filter and return the filter coefficients
    ///
    /// # Parameters
    ///  - order : Order of the filter
    ///  - band: [`BandFilter`] to apply, for analog filters band is expressed as an angular</br>
    /// frequency (rad/s).
    /// for digital filters band is in the same units as [`Analog`]
    ///  - analog: Analog or Digital filter selection, when digital is selected<br/>
    ///  a sampling rate is required
    ///  - norm: Critical frequency normalization, see [`BesselNorm`] for more info
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog, norm: BesselNorm) -> Ba {
        bessel_filter(order, band_filter, analog, norm, DesiredFilterOutput::Ba).ba()
    }
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

/// TODO! _norm defaults to Phase, other normalizations are not implemented
pub fn besselap(order: u32, _norm: BesselNorm) -> Zpk {
    let z = array![];
    let mut p: Array1<Complex<f64>>;
    let k = 1.0;
    if order == 0 {
        p = array![];
    } else {
        p = _bessel_zeros(order).into_iter().map(|a| 1.0 / a).collect();
        let a_last = (_falling_factorial(2 * order, order) / 2.0.powi(order as _)).floor();
        p.iter_mut()
            .for_each(|a| *a *= 10.0.powf(-a_last.log10() / order as f64));
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

fn _bessel_zeros(order: u32) -> Array1<Complex<f64>> {
    if order == 0 {
        return Default::default();
    }

    let x0 = _campos_zeros(order);
    let f = |x: Complex<f64>| kve(order as f64 + 0.5, 1.0 / x);

    let fp = |x: Complex<f64>| {
        let order = order as f64;

        let first = kve(order - 0.5, 1.0 / x) / (2.0 * x.powi(2));
        let second = kve(order + 0.5, 1.0 / x) / x.powi(2);
        let third = kve(order + 1.5, 1.0 / x) / (2.0 * x.powi(2));
        first - second + third
    };
    let mut x = _aberth(f, fp, &x0);

    for i in &mut x {
        *i = newton(f, fp, *i, 10.0.powi(-16), 50);
    }

    let clone = x.clone().into_iter().map(|a| a.conj()).rev();

    let temp = x.iter().copied().zip(clone);
    let x: Array1<Complex<f64>> = temp.map(|(a, b)| (a + b) / 2.0).collect();

    x
}

fn _aberth<F: Fn(Complex<f64>) -> Complex<f64>, FP: Fn(Complex<f64>) -> Complex<f64>>(
    f: F,
    fp: FP,
    x0: &[Complex<f64>],
) -> Vec<Complex<f64>> {
    let mut zs = x0.to_vec();
    let mut new_zs = zs.clone();
    let tol = 10.0.powi(-16);
    'iteration: for _ in 0..100 {
        for i in 0..(x0.len()) {
            let p_of_z = f(zs[i]);
            let dydx_of_z = fp(zs[i]);

            let sum: Complex64 = (0..zs.len())
                .filter(|&k| k != i)
                .fold(Complex::zero(), |acc: Complex64, k| {
                    acc + Complex64::one() / (zs[i] - zs[k])
                });

            let new_z = zs[i] + p_of_z / (p_of_z * sum - dydx_of_z);
            new_zs[i] = new_z;
            if new_z.re.is_nan()
                || new_z.im.is_nan()
                || new_z.re.is_infinite()
                || new_z.im.is_infinite()
            {
                break 'iteration;
            }
            let err = (new_z - zs[i]).abs();
            if err < tol {
                return new_zs;
            }

            zs = new_zs.clone();
        }
    }

    panic!();
}

// verified with python
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

    assert_eq!(x.len(), y.len());
    x.iter()
        .zip(y)
        .map(|(x, y)| *x + Complex::new(0.0, 1.0) * y)
        .collect::<Vec<_>>()
}

#[cfg(test)]
mod tests {
    use num::Complex;

    use crate::{
        signal::{
            band_filter::BandFilter,
            bessel::{_aberth, _campos_zeros},
            output_type::Zpk,
            tools::newton,
        },
        special::kve,
    };

    use super::{besselap, BesselFilterStandalone, BesselNorm};

    #[test]
    fn test_bessel_filter() {
        let _filter = BesselFilterStandalone::<Zpk>::filter(
            4,
            BandFilter::Lowpass(0.2),
            crate::signal::Analog::False { fs: 2.0 },
            BesselNorm::Phase,
        );
    }

    #[test]
    fn test_besselap() {
        let _res = besselap(4, BesselNorm::Phase);
    }
    #[test]
    fn test_besselzeros() {
        let res = _campos_zeros(4);
        let order = 4;
        let f = |x: Complex<f64>| kve(order as f64 + 0.5, 1.0 / x);

        let fp = |x: Complex<f64>| {
            let order = order as f64;

            let first = kve(order - 0.5, 1.0 / x) / (2.0 * x.powi(2));
            let second = kve(order + 0.5, 1.0 / x) / x.powi(2);
            let third = kve(order + 1.5, 1.0 / x) / (2.0 * x.powi(2));
            first - second + third
        };
        let mut res = _aberth(f, fp, &res);
        for i in &mut res {
            *i = newton(f, fp, *i, 10.0_f64.powi(-16), 50);
        }
    }
}
