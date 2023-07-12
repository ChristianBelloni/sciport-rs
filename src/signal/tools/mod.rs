use std::ops::{Div, Mul, Sub};

use num::{
    complex::{Complex64, ComplexFloat},
    Complex, Num, One, Signed, Zero,
};

use super::output_type::Zpk;

pub fn bilinear_zpk(mut input: Zpk, fs: f64) -> Zpk {
    println!("z size {} p size {}", input.z.len(), input.p.len());
    let degree = relative_degree(&input);
    let fs2 = fs * 2.0;

    println!("fs - {fs}");
    let z_prod = input
        .z
        .iter()
        .fold(Complex::one(), |acc: Complex<f64>, item| acc * (fs2 - item));
    let p_prod = input
        .p
        .iter()
        .fold(Complex::one(), |acc: Complex<f64>, item| acc * (fs2 - item));

    input.z = input.z.iter().map(|&a| (fs2 + a) / (fs2 - a)).collect();
    input.p = input.p.iter().map(|&a| (fs2 + a) / (fs2 - a)).collect();

    input.z.extend(vec![-Complex::one(); degree]);
    println!("{z_prod} -- {p_prod}");
    input.k *= (z_prod / p_prod).re;
    input
}

pub fn relative_degree(input: &Zpk) -> usize {
    input.p.len() - input.z.len()
}

pub fn poly<T: Num + Copy>(zeroes: &[T]) -> Vec<T> {
    let mut coeff = vec![T::one()];
    for z in zeroes {
        let mut clone = coeff.clone();
        mul_by_x(&mut coeff);
        mul_by_scalar(&mut clone, *z);
        sub_coeff(&mut coeff[1..], &clone)
    }
    coeff
}

fn mul_by_x<T: Num>(coeff: &mut Vec<T>) {
    coeff.push(T::zero());
}

fn mul_by_scalar<T: Num + Copy>(coeff: &mut [T], scalar: T) {
    coeff.iter_mut().for_each(move |a| *a = *a * scalar);
}

fn sub_coeff<T: Num + Copy>(coeff: &mut [T], ar: &[T]) {
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

pub fn find_root<N>(
    function: impl Fn(N) -> N,
    derivative: impl Fn(N) -> N,
    x0: N,
    _acceptable_err: N,
    max_iterations: i32,
) -> Result<N, N>
where
    N: Div<Output = N> + Sub<Output = N>,
    N: One<Output = N> + Copy + Mul<N, Output = N>,
{
    let mut current_x: N = x0;
    let mut next_x: N;

    for _ in 0..max_iterations {
        let deviation = function(current_x) / derivative(current_x);
        next_x = current_x - deviation;
        current_x = next_x;
    }
    Ok(current_x)
}
pub fn find_root_complex(
    f: impl Fn(Complex64) -> Complex64,
    fp: impl Fn(Complex64) -> Complex64,
    x0: Complex64,
    acceptable_err: f64,
    max_iterations: i32,
) -> Complex64 {
    let mut current_x = x0;
    let mut next_x;

    for _ in 0..max_iterations {
        let deviation = -f(current_x) / fp(current_x);

        next_x = current_x + deviation;
        current_x = next_x;
        if deviation.abs() < acceptable_err {
            return current_x;
        }
    }
    panic!();
}

pub fn newton(
    f: impl Fn(Complex64) -> Complex64,
    fp: impl Fn(Complex64) -> Complex64,
    mut x0: Complex64,
    tol: f64,
    maxiter: usize,
) -> Complex64 {
    for _ in 0..maxiter {
        let fval = f(x0);

        if fval == 0.0.into() {}

        let fder = fp(x0);
        let newton_step = fval / fder;
        let x = x0 - newton_step;

        if is_close(x, x0, tol) {
            return x;
        }
        x0 = x;
    }

    panic!()
}

fn is_close(x: Complex64, y: Complex64, tol: f64) -> bool {
    let df = x - y;
    let a = df.abs();

    a < tol
}

#[cfg(test)]
mod tests {
    use num::complex::Complex64;

    use crate::special::kv::kve;
    #[test]
    fn test_i() {
        let res = kve(1.0, Complex64::new(-29.5, -88.33333333));
        dbg!(res);
    }
}
