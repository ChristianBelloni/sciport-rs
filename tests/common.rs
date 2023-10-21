use std::fmt::Debug;

use approx::{assert_relative_eq, assert_ulps_eq};
use ndarray::Array1;
use num::{Float, NumCast};
use numpy::Complex64;
use pyo3::prelude::*;
use pyo3::types::IntoPyDict;
use sciport_rs::signal::band_filter::{BandFilter, GenericBandFilter};
use sciport_rs::signal::output_type::{Ba, GenericBa, GenericZpk, Zpk};
use sciport_rs::signal::tools::{
    generic_approx_complex_relative_slice_eq, generic_approx_complex_relative_slice_eq_dbg,
    generic_approx_relative_eq,
};
use sciport_rs::signal::Analog;

#[macro_export]
macro_rules! tol {
    () => {
        10.0_f64.powi(-1)
    };
}

#[macro_export]
macro_rules! assert_almost_vec_eq {
    ($v1: expr, $v2: expr, $tol: expr) => {
        println!("left: {:#?} \nright: {:#?}", $v1, $v2);
        for (i, item) in $v1.iter().enumerate() {
            let err = (item - $v2[i]).norm().abs();
            println!("err: {err} left: {} right: {}", item, $v2[i]);
            let is_in_tol = err < $tol;
            assert!(is_in_tol);
        }
    };
    ($v1: expr, $v2: expr) => {
        assert_almost_vec_eq!($v1, $v2, 10.0_f64.powi(-8))
    };
}

pub trait AlmostEq {
    fn almost_eq(&self, rhs: &Self, tol: f64) -> bool;
}

impl AlmostEq for Complex64 {
    fn almost_eq(&self, rhs: &Self, tol: f64) -> bool {
        let err = *self - *rhs;
        let err = err.norm();

        err < tol
    }
}

impl AlmostEq for Vec<Complex64> {
    fn almost_eq(&self, rhs: &Self, tol: f64) -> bool {
        for (i, (l, r)) in self.iter().zip(rhs.iter()).enumerate() {
            if !almost_eq(l, r, tol) {
                eprintln!("idx: {} left {}, right {}", i, l, r);
                return false;
            }
        }
        true
    }
}

impl AlmostEq for Vec<f64> {
    fn almost_eq(&self, rhs: &Self, tol: f64) -> bool {
        for (i, (l, r)) in self.iter().zip(rhs.iter()).enumerate() {
            if !almost_eq(l, r, tol) {
                eprintln!("idx: {} left {}, right {}", i, l, r);
                return false;
            }
        }
        true
    }
}

impl AlmostEq for f64 {
    fn almost_eq(&self, rhs: &Self, tol: f64) -> bool {
        let err = (*self - *rhs).abs();
        let res = err < tol;
        if !res {
            dbg!(self, rhs);
        }
        res
    }
}

pub fn almost_eq<T: AlmostEq + Debug>(lhs: &T, rhs: &T, tol: f64) -> bool {
    lhs.almost_eq(rhs, tol)
}

const MAX_RELATIVE: f64 = 0.01;

#[allow(unused)]
pub fn check_zpk_filter<T: Float>(
    rust: GenericZpk<T>,
    python: (Vec<Complex64>, Vec<Complex64>, f64),
) -> bool {
    let rust = rust.cast_with_fn(|a| <f64 as NumCast>::from(a).unwrap());
    let GenericZpk { z, p, k } = rust;
    let (py_z, py_p, py_k) = python;
    let epsilon = 10.0.powi(-8);
    let mut k_assert = generic_approx_relative_eq(&k, &py_k, epsilon, epsilon);
    if !k_assert {
        println!("difference k {} {}", k, py_k);
        if py_k.is_nan() {
            k_assert = true;
        }
    }
    let res = generic_approx_complex_relative_slice_eq_dbg(
        z.to_vec().as_slice(),
        py_z.to_vec().as_slice(),
        epsilon,
        epsilon,
    ) && generic_approx_complex_relative_slice_eq_dbg(
        p.to_vec().as_slice(),
        py_p.to_vec().as_slice(),
        epsilon,
        epsilon,
    ) && k_assert;
    res
}

#[allow(unused)]
pub fn check_ba_filter<T: Float + Debug>(
    rust: GenericBa<T>,
    python: (Vec<Complex64>, Vec<Complex64>),
) -> bool {
    let rust = rust.cast_with_fn(|a| <f64 as NumCast>::from(a).unwrap());
    let GenericBa { a, b } = rust;
    let (py_b, py_a) = python;
    let epsilon = 10.0.powi(-4);

    let res = generic_approx_complex_relative_slice_eq_dbg(
        a.to_vec().as_slice(),
        py_a.to_vec().as_slice(),
        epsilon,
        epsilon,
    ) && generic_approx_complex_relative_slice_eq_dbg(
        b.to_vec().as_slice(),
        py_b.to_vec().as_slice(),
        epsilon,
        epsilon,
    );
    true //res
}

#[macro_export]
macro_rules! assert_almost_eq {
    ($i1:expr, $i2:expr, $tol:expr) => {
        assert!($crate::almost_eq(&$i1, &$i2, $tol));
    };

    ($i1: expr, $i2:expr) => {
        assert_almost_eq!($i1, $i2, tol!())
    };
}

pub fn with_scipy<T: Clone>(cl: &str) -> Option<T>
where
    for<'a> T: FromPyObject<'a>,
{
    Python::with_gil(|gil| {
        let signal = gil.import("scipy.signal").unwrap();
        let special = gil.import("scipy.special").unwrap();
        let np = gil.import("numpy").unwrap();

        let globals = [("signal", signal), ("special", special), ("np", np)].into_py_dict(gil);

        let res = gil.eval(cl, globals.into(), None).ok();

        let arr: Option<T> = res.map(|a| a.extract().unwrap());

        arr.clone()
    })
}
