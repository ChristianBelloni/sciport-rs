use std::ffi::{c_double, c_int};

extern "C" {

    fn zbesy_wrap(
        zr: c_double,
        zi: c_double,
        nu: c_double,
        kode: c_int,
        N: c_int,
        cyr: *mut c_double,
        cyi: *mut c_double,
        nz: *mut c_int,
        cwrkr: *mut c_double,
        cwrki: *mut c_double,
        ierr: *mut c_int,
    );

    fn zbesj_wrap(
        zr: c_double,
        zi: c_double,
        nu: c_double,
        kode: c_int,
        n: c_int,
        cyr: *mut c_double,
        cyi: *mut c_double,
        nz: *mut c_int,
        ierr: *mut c_int,
    );
}

pub use wrappers::*;

mod wrappers {
    use num::complex::Complex64;

    use super::{zbesj_wrap, zbesy_wrap};

    pub fn bessel_y(order: f64, z: Complex64) -> Result<Complex64, i32> {
        unsafe { _bessel_y(order, z) }
    }

    unsafe fn _bessel_y(order: f64, z: Complex64) -> Result<Complex64, i32> {
        let zr = z.re;
        let zi = z.im;
        let nu = order.abs();
        let kode = 1;
        let n = 1;

        let mut cyr = 0.0;
        let mut cyi = 0.0;
        let mut cwrkr = 0.0;
        let mut cwrki = 0.0;
        
        let mut nz = 0;
        let mut ierr = 0;

        let mut answer;

        zbesy_wrap(
            zr, zi, nu, kode, n, &mut cyr, &mut cyi, &mut nz, &mut cwrkr, &mut cwrki, &mut ierr,
        );

        if zi == 0.0 && zr == 0.0 {
            cyi = 0.0;
        }

        answer = Complex64::new(cyr, cyi);

        if ierr != 0 {
            Err(ierr)?;
        }

        if order < 0.0 {
            let c = nu.cos();
            let s = nu.sin();

            let mut cyrj = 0.0;
            let mut cyij = 0.0;

            let mut nz_j = 0;
            let mut ierrj = 0;

            let kodej = 1;
            let nj = 1;

            zbesj_wrap(
                zr, zi, nu, kodej, nj, &mut cyrj, &mut cyij, &mut nz_j, &mut ierrj,
            );

            if ierrj != 0 {
                Err(ierrj)?;
            }

            let answer_j = Complex64::new(cyrj, cyij);

            answer = s * answer_j + c * answer;

            return Ok(answer);
        }

        Ok(answer)
    }
}


#[cfg(test)]
mod tests {
    use super::wrappers::bessel_y;
    use num::Complex;

    #[test]
    fn test_ffi() {
        let answer = bessel_y(-4.0, Complex::new(2.0, 1.0));

        println!("{answer:#?}")
    }
}
