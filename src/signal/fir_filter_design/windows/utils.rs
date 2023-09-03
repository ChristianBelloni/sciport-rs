use ndarray::Array1;

pub fn len_guards(m: u64) -> bool {
    m <= 1
}

pub fn extend(m: u64, sym: bool) -> (u64, bool) {
    if !sym {
        (m + 1, true)
    } else {
        (m, false)
    }
}

pub fn truncate(w: Array1<f64>, needed: bool) -> Array1<f64> {
    if needed {
        let mut inner = w.to_vec();
        inner.truncate(inner.len() - 1);
        Array1::from_shape_vec(inner.len(), inner).unwrap()
    } else {
        w
    }
}
