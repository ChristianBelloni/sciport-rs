use num::Complex;
use pyo3::{types::IntoPyDict, FromPyObject, Python};
use rand::Rng;

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
#[test]
pub fn test_i0() {
    for _ in 0..500 {
        let len = rand::thread_rng().gen_range(2..100);
        let mut v = Vec::with_capacity(len);
        for _ in 0..len {
            v.push(rand::thread_rng().gen_range(0.0..1.0));
        }

        let rust_result = sciport_rs::special::i0(v.clone().into())
            .unwrap()
            .mapv(|a| a.norm())
            .to_vec();
        let py_script = format!("special.i0({v:?})");
        let py_res: Vec<Complex<f64>> = with_scipy(&py_script).unwrap();
        let py_res: Vec<f64> = py_res.into_iter().map(|a| a.norm()).collect();

        approx::assert_relative_eq!(
            rust_result.as_slice(),
            py_res.as_slice(),
            epsilon = 0.00000000000001
        );
    }
}
