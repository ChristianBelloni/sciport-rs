use num::{complex::Complex64, Complex, Num, One, Zero};

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
    input.k = input.k * (z_prod / p_prod).re;
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

fn mul_by_scalar<T: Num + Copy>(coeff: &mut Vec<T>, scalar: T) {
    coeff.iter_mut().for_each(move |a| *a = *a * scalar);
}

fn sub_scalar<T: Num + Copy>(coeff: &mut Vec<T>, scalar: T) {
    coeff.iter_mut().for_each(|a| *a = *a - scalar);
}

fn sub_coeff<T: Num + Copy>(coeff: &mut [T], ar: &[T]) {
    for (i, c) in coeff.iter_mut().enumerate() {
        *c = *c - ar[i]
    }
}

pub fn polyval<T: Into<Complex<f64>> + Copy, const S: usize>(v: T, coeff: [T; S]) -> Complex<f64> {
    fn polyval<const S: usize>(v: Complex<f64>, coeff: [Complex<f64>; S]) -> Complex<f64> {
        let max_exp = coeff.len();

        let it = coeff
            .iter()
            .enumerate()
            .map(|(i, c)| {
                let exp = max_exp - i;

                v.powi(exp as _) * c
            })
            .sum();
        it
    }

    let tmp: [Complex<f64>; S] = coeff
        .into_iter()
        .map(|a| a.into())
        .collect::<Vec<_>>()
        .try_into()
        .unwrap();
    polyval(v.into(), tmp.into())
}


#[cfg(test)]
mod tests {
    use num::complex::Complex64;

    use crate::special::kv::{kv, kve};
    #[test]
    fn test_i() {
        let res = kve(1.0, Complex64::new(-29.5, -88.33333333));
        dbg!(res);
    }
}
