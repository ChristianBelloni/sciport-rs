use std::marker::PhantomData;

use super::output_type::Zpk;
use num::complex::Complex64;
use num::One;

pub struct Cheb1FilterStandalone<T>(PhantomData<T>);

impl Cheb1FilterStandalone<Zpk> {
    pub fn filter(order: u32, band: BandFilter)
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
        k = k / (1.0 + eps * eps).sqrt();
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
