use super::output_type::{GenericBa, GenericZpk};
use crate::odr::polynomial::Polynomial;
use ndarray::{array, concatenate, s, Array1, ArrayView1, ArrayViewMut1, Axis};
use num::{
    complex::{Complex64, ComplexFloat},
    traits::FloatConst,
    Complex, Float, Num, Zero,
};
use std::fmt::Debug;
use trait_set::trait_set;

trait_set! {
    pub trait BilinearZpk = Float;
}

pub fn bilinear_zpk<T>(input: GenericZpk<T>, fs: T) -> GenericZpk<T>
where
    T: BilinearZpk,
{
    let degree = relative_degree(&input);
    let GenericZpk { z, p, k } = input;

    let from = |a: f64| T::from(a).unwrap();

    let fs2 = from(2.0) * fs;

    let z_z = z.mapv(|a| a + fs2) / z.mapv(|a| -a + fs2);
    let p_z = p.mapv(|a| a + fs2) / p.mapv(|a| -a + fs2);

    let z_z = concatenate![Axis(0), z_z, -Array1::ones(degree)];

    let factor_map = match z.len().cmp(&p.len()) {
        std::cmp::Ordering::Equal => z.mapv(|a| -a + fs2) / p.mapv(|a| -a + fs2),
        std::cmp::Ordering::Greater => {
            let t_z = z.mapv(|a| -a + fs2);
            let t_p = p.mapv(|a| -a + fs2);
            let t_p = concatenate![Axis(0), t_p, Array1::ones(z.len() - p.len())];
            t_z / t_p
        }
        std::cmp::Ordering::Less => {
            let t_z = z.mapv(|a| -a + fs2);
            let t_p = p.mapv(|a| -a + fs2);
            let t_z = concatenate![Axis(0), t_z, Array1::ones(p.len() - z.len())];
            t_z / t_p
        }
    };

    let k_z = k * factor_map.product().re;

    GenericZpk {
        z: z_z,
        p: p_z,
        k: k_z,
    }
}

pub fn relative_degree<T>(input: &GenericZpk<T>) -> usize {
    input.p.len() - input.z.len()
}

/// Compute polynomial coefficients from zeroes
///
/// # Examples
///
/// ```rust
/// # use sciport_rs::signal::tools::poly;
/// # use num::complex::Complex64;
/// # use ndarray::array;
/// let complex_z = array![ Complex64::new(2.1, 3.2), Complex64::new(1.0, 1.0) ];
///
/// let coeffs = poly((&complex_z).into());
///
/// ```
pub fn poly<T: Num + Copy>(zeroes: ArrayView1<'_, T>) -> Array1<T> {
    let mut coeff = array![T::one()];
    for z in zeroes {
        let mut clone = coeff.clone();
        mul_by_x(&mut coeff);
        mul_by_scalar(clone.view_mut(), *z);
        sub_coeff(coeff.slice_mut(s![1..]), clone)
    }
    coeff
}

fn mul_by_x<T: Num + Clone>(coeff: &mut Array1<T>) {
    coeff.append(Axis(0), (&array![T::zero()]).into()).unwrap();
}

fn mul_by_scalar<T: Num + Copy>(mut coeff: ArrayViewMut1<T>, scalar: T) {
    coeff.map_inplace(move |a| *a = *a * scalar);
}

fn sub_coeff<T: Num + Copy>(mut coeff: ArrayViewMut1<T>, ar: Array1<T>) {
    for (i, c) in coeff.iter_mut().enumerate() {
        *c = *c - ar[i]
    }
}

pub fn polyval<T: Into<Complex<f64>> + Copy, const S: usize>(v: T, coeff: [T; S]) -> Complex<f64> {
    fn polyval<const S: usize>(v: Complex<f64>, coeff: [Complex<f64>; S]) -> Complex<f64> {
        coeff
            .iter()
            .enumerate()
            .fold(Complex64::zero(), |acc, (i, item)| {
                acc + v.powi(i as _) * item
            })
    }

    let tmp: [Complex<f64>; S] = coeff
        .into_iter()
        .map(|a| a.into())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    polyval(v.into(), tmp)
}

pub fn zpk2ba<T>(zpk: GenericZpk<T>) -> GenericBa<T>
where
    T: Float + FloatConst + ComplexFloat,
{
    let GenericZpk { z, p, k } = zpk;

    let pol = Polynomial::from_roots_k(z.clone(), T::one().into());
    let pol = pol.saturate();
    let mut b: Array1<_> = pol.iter().rev().copied().collect();

    b = b.mapv(|a| a * k);

    let mut a: Array1<_> = Polynomial::from_roots_k(p.clone(), T::one().into())
        .saturate()
        .into_iter()
        .rev()
        .collect();
    let epsilon = T::from(10.0_f64.powi(-4)).unwrap();

    let roots = z;
    let mut pos_roots: Vec<Complex<T>> =
        roots.iter().cloned().filter(|a| a.im > T::zero()).collect();
    let mut neg_roots: Vec<Complex<T>> = roots
        .iter()
        .filter(|a| a.im < T::zero())
        .map(|a| a.conj())
        .collect();

    if pos_roots.len() == neg_roots.len() {
        sort_complex(&mut pos_roots);
        sort_complex(&mut neg_roots);

        if generic_approx_complex_relative_slice_eq(
            pos_roots.as_slice(),
            neg_roots.as_slice(),
            epsilon,
            epsilon,
        ) {
            b = b.into_iter().map(|a| a.re.into()).collect();
        }
    }

    let roots = p;
    let mut pos_roots: Vec<Complex<T>> =
        roots.iter().cloned().filter(|a| a.im > T::zero()).collect();
    let mut neg_roots: Vec<Complex<T>> = roots
        .iter()
        .filter(|a| a.im < T::zero())
        .map(|a| a.conj())
        .collect();

    if pos_roots.len() == neg_roots.len() {
        sort_complex(&mut pos_roots);
        sort_complex(&mut neg_roots);

        // println!("pos_roots: {:?}", pos_roots);
        // println!("neg_roots: {:?}", neg_roots);

        if generic_approx_complex_relative_slice_eq(
            pos_roots.as_slice(),
            neg_roots.as_slice(),
            epsilon,
            epsilon,
        ) {
            a = a.into_iter().map(|a| a.re.into()).collect();
        }
    }
    // this shouldn't be here but without it nothing works!
    if b[0] == Complex::zero() {
        b.remove_index(Axis(0), 0);
    }

    if a[0] == Complex::zero() {
        a.remove_index(Axis(0), 0);
    }
    let a = a;
    let b = b;
    GenericBa { a, b }
}

/// lifted from https://docs.rs/approx/latest/src/approx/relative_eq.rs.html#8-30
pub fn generic_approx_relative_eq<T: Float + Clone, K: Float + Clone>(
    lhs: &T,
    rhs: &K,
    epsilon: T,
    max_relative: T,
) -> bool {
    let rhs = *rhs;
    let lhs = *lhs;
    if lhs == T::from(rhs).unwrap() {
        return true;
    }
    if lhs.is_infinite() || rhs.is_infinite() {
        return false;
    }

    let abs_diff = (lhs - T::from(rhs).unwrap()).abs();
    if abs_diff < epsilon {
        return true;
    }

    let abs_lhs = lhs.abs();
    let abs_rhs = rhs.abs();

    let largest = if abs_lhs > T::from(abs_rhs).unwrap() {
        abs_lhs
    } else {
        T::from(abs_rhs).unwrap()
    };

    abs_diff <= largest * max_relative
}
pub fn generic_approx_complex_relative_eq<T: Float + Clone, K: Float + Clone>(
    lhs: &Complex<T>,
    rhs: &Complex<K>,
    epsilon: T,
    max_relative: T,
) -> bool {
    generic_approx_relative_eq(&lhs.re, &rhs.re, epsilon, max_relative)
        && generic_approx_relative_eq(&lhs.im, &rhs.im, epsilon, max_relative)
}

pub fn generic_approx_complex_relative_eq_dbg<
    T: Float + Clone + Debug,
    K: Float + Clone + Debug,
>(
    lhs: &Complex<T>,
    rhs: &Complex<K>,
    epsilon: T,
    max_relative: T,
) -> bool {
    let res = generic_approx_relative_eq(&lhs.re, &rhs.re, epsilon, max_relative)
        && generic_approx_relative_eq(&lhs.im, &rhs.im, epsilon, max_relative);
    if !res {
        println!("difference {:?} {:?}", lhs, rhs);
    }
    res
}

pub fn generic_approx_relative_slice_eq<T: Float + Clone>(
    lhs: &[T],
    rhs: &[T],
    epsilon: T,
    max_relative: T,
) -> bool {
    let zip = lhs.iter().zip(rhs.iter());
    zip.fold(true, |acc, (lhs, rhs)| {
        acc && generic_approx_relative_eq(lhs, rhs, epsilon, max_relative)
    })
}

pub fn generic_approx_complex_relative_slice_eq<T: Float + Clone, K: Float + Clone>(
    lhs: &[Complex<T>],
    rhs: &[Complex<K>],
    epsilon: T,
    max_relative: T,
) -> bool {
    let zip = lhs.iter().zip(rhs.iter());
    zip.enumerate().fold(true, |acc, (i, (lhs, rhs))| {
        let new = generic_approx_complex_relative_eq(lhs, rhs, epsilon, max_relative);
        if !new {
            println!("difference at {}", i);
        }
        acc && new
    })
}

pub fn generic_approx_complex_relative_slice_eq_dbg<
    T: Float + Clone + Debug,
    K: Float + Clone + Debug,
>(
    lhs: &[Complex<T>],
    rhs: &[Complex<K>],
    epsilon: T,
    max_relative: T,
) -> bool {
    let zip = lhs.iter().zip(rhs.iter());
    zip.fold(true, |acc, (lhs, rhs)| {
        let new = generic_approx_complex_relative_eq_dbg(lhs, rhs, epsilon, max_relative);

        acc && new
    })
}

fn sort_complex<T>(cxs: &mut [Complex<T>])
where
    T: Float,
{
    cxs.sort_by(|a, b| match a.re.partial_cmp(&b.re) {
        Some(comp) => match comp {
            std::cmp::Ordering::Less => std::cmp::Ordering::Less,
            std::cmp::Ordering::Greater => std::cmp::Ordering::Greater,
            std::cmp::Ordering::Equal => {
                a.im.partial_cmp(&b.im).unwrap_or(std::cmp::Ordering::Equal)
            }
        },
        None => a.im.partial_cmp(&b.im).unwrap_or(std::cmp::Ordering::Equal),
    })
}
