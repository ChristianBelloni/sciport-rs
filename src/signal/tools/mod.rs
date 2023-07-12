use num::{
    complex::{Complex64, ComplexFloat},
    Complex, Num, One, Zero,
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

/// Compute polynomial coefficients from zeroes
///
/// # Examples
///
/// ```rust
/// # use sciport_rs::signal::tools::poly;
/// # use num::complex::Complex64;
/// let zeros = [2, 3];
///
/// let coeffs = poly(&zeros);
///
/// assert_eq!(&coeffs, &[1, -5, 6]);
///
/// let complex_z = [ Complex64::new(2.1, 3.2), Complex64::new(1.0, 1.0) ];
///
/// let coeffs = poly(&complex_z);
///
/// ```
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

    use crate::special::kve;
    #[test]
    fn test_i() {
        let res = kve(1.0, Complex64::new(-29.5, -88.33333333));
        dbg!(res);
    }
}
