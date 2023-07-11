use num::{complex::Complex64, Complex, One};
pub fn kv(v: f64, mut z: Complex<f64>) -> Complex<f64> {
    if z.is_nan() {
        panic!();
    }
    complex_bessel_rs::bessel_k::bessel_k(v, z).unwrap()
}

pub fn kve(v: f64, z: Complex<f64>) -> Complex<f64> {
    kv(v, z) * z.exp()
}

#[cfg(test)]
mod tests {
    use num::complex::Complex64;

    use crate::special::kv::kve;

    use super::kv;

    #[test]
    fn test_kv() {
        let c1 = Complex64::new(1000.0, 0.0);
        let res = kv(-3.5, c1);
        let res2 = kve(0.0, 1.0.into());
        let res3 = kve(3.5, 1.0 / Complex64::new(1.2, 0.3));
        println!("{res}");
        println!("{res2}");
        println!("{res3}");
    }
}
