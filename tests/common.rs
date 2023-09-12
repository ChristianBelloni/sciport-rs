use std::fmt::Debug;

use numpy::Complex64;
use pyo3::prelude::*;
use pyo3::types::IntoPyDict;
use sciport_rs::signal::output_type::Zpk;

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

pub fn check_zpk_filter(rust: Zpk, python: (Vec<Complex64>, Vec<Complex64>, f64)) -> bool {
    let Zpk { z, p, k } = rust;
    let (py_z, py_p, py_k) = python;

    let res = almost_eq(&py_z, &z.to_vec(), tol!())
        && almost_eq(&py_p, &p.to_vec(), tol!())
        && almost_eq(&py_k, &k, tol!());
    if !res && false {
        println!("rust {z:#?} {p:#?} {k:#?}");
        println!("python {py_z:#?} {py_p:#?} {py_k:#?}");
    }
    res
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

pub fn with_scipy<T: Clone>(cl: &str) -> T
where
    for<'a> T: FromPyObject<'a>,
{
    Python::with_gil(|gil| {
        let signal = gil.import("scipy.signal").unwrap();

        let globals = [("signal", signal)].into_py_dict(gil);

        let res = gil.eval(cl, globals.into(), None).unwrap();

        let arr: T = res.extract().unwrap();

        arr.clone()
    })
}
