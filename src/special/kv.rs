use num::Complex;
pub fn kv(v: f64, z: Complex<f64>) -> Complex<f64> {
    crate::bindings::bessel_y(v, z).unwrap()
}

pub fn kve(v: f64, z: Complex<f64>) -> Complex<f64> {
    (kv(v, z)) * (z.exp())
}

#[cfg(test)]
mod tests {
    use num::complex::Complex64;

    use super::kv;

    #[test]
    fn test_kv() {
        let c1 = Complex64::new(2.0, -1.0);
        let res = kv(-3.5, c1);
        println!("{res}");
    }
}
