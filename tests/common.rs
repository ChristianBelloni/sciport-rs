use pyo3::prelude::*;

pub fn with_scipy(cl: impl Fn(Python<'_>)) {
    Python::with_gil(cl);
}
