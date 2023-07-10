use std::cmp::Ordering;

use num::{Complex, Zero};

use crate::signal::tools::poly;

use super::{Ba, Zpk};

impl From<Zpk> for Ba {
    fn from(value: Zpk) -> Self {
        let mut b = poly(&value.z);

        b.iter_mut().for_each(|a| *a *= value.k);

        let mut a = poly(&value.p);
        {
            let mut pos_z_roots = value
                .z
                .iter()
                .filter(|a| a.im > 0.0)
                .copied()
                .collect::<Vec<_>>();
            let mut neg_z_roots = value
                .z
                .iter()
                .filter(|a| a.im < 0.0)
                .map(|a| a.conj())
                .collect::<Vec<_>>();

            if pos_z_roots.len() == neg_z_roots.len() {
                pos_z_roots.sort_by(|a, b| match a.re.total_cmp(&b.re) {
                    Ordering::Less => Ordering::Less,
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Equal => a.im.total_cmp(&b.im),
                });

                neg_z_roots.sort_by(|a, b| match a.re.total_cmp(&b.re) {
                    Ordering::Less => Ordering::Less,
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Equal => a.im.total_cmp(&b.im),
                });

                if pos_z_roots == neg_z_roots {
                    b.iter_mut().for_each(|a| *a = a.re.into());
                }
            }
        }
        {
            let mut pos_p_roots = value
                .p
                .iter()
                .filter(|a| a.im > 0.0)
                .copied()
                .collect::<Vec<_>>();
            let mut neg_p_roots = value
                .p
                .iter()
                .filter(|a| a.im < 0.0)
                .map(|a| a.conj())
                .collect::<Vec<_>>();

            if pos_p_roots.len() == neg_p_roots.len() {
                pos_p_roots.sort_by(|a, b| match a.re.total_cmp(&b.re) {
                    Ordering::Less => Ordering::Less,
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Equal => a.im.total_cmp(&b.im),
                });

                neg_p_roots.sort_by(|a, b| match a.re.total_cmp(&b.re) {
                    Ordering::Less => Ordering::Less,
                    Ordering::Greater => Ordering::Greater,
                    Ordering::Equal => a.im.total_cmp(&b.im),
                });

                if pos_p_roots == neg_p_roots {
                    a.iter_mut().for_each(|a| *a = a.re.into());
                }
            }
        }
        Self { a, b }
    }
}
#[allow(unused)]
fn mul_by_x(input: &mut Vec<Complex<f64>>) {
    input.push(Complex::zero());
}

#[allow(unused)]
fn mul_by_scalar(input: &mut [Complex<f64>], scalar: Complex<f64>) {
    input.iter_mut().for_each(|a| *a *= scalar);
}

#[allow(unused)]
fn sub(input: &mut [Complex<f64>], other: &[Complex<f64>]) {
    for (i, item) in other.iter().enumerate() {
        *input.get_mut(i).unwrap() = *input.get(i).unwrap() - item;
    }
}
