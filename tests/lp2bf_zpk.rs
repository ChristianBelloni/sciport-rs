mod common;

use crate::common::check_zpk_filter;
use common::with_scipy;
use num::complex::Complex64;
use rand::Rng;
use sciport_rs::signal::{
    band_filter::{lp2bp_zpk, lp2bs_zpk, lp2hp_zpk, lp2lp_zpk},
    butter::*,
};

#[test]
fn with_py_fuzzy_lp2lp_zpk() {
    for _ in 0..10_000 {
        let mut rng = rand::thread_rng();
        let order = rng.gen_range(1..200);
        let input = buttap(order as _);
        let wo = rng.gen_range(0.0..1.0);
        let input2 = input.clone();
        let rust = lp2lp_zpk(input.clone(), wo);

        let z_formatted = format!("{}", input2.z).replace("i", "j");
        let p_formatted = format!("{}", input2.p).replace("i", "j");

        let py_code = format!(
            "signal.lp2lp_zpk({z_formatted},{p_formatted},{}, wo={wo})",
            input.k
        );

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);

        assert!(check_zpk_filter(rust, python))
    }
}

#[test]
fn with_py_fuzzy_lp2hp_zpk() {
    for _ in 0..10_000 {
        let mut rng = rand::thread_rng();
        let order = rng.gen_range(1..200);
        let input = buttap(order as _);
        let wo = rng.gen_range(0.0..1.0);
        let input2 = input.clone();
        let rust = lp2hp_zpk(input.clone(), wo);

        let z_formatted = format!("{}", input2.z).replace("i", "j");
        let p_formatted = format!("{}", input2.p).replace("i", "j");

        let py_code = format!(
            "signal.lp2hp_zpk({z_formatted},{p_formatted},{}, wo={wo})",
            input.k
        );

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);

        assert!(check_zpk_filter(rust, python))
    }
}

#[test]
fn with_py_fuzzy_lp2bp_zpk() {
    for _ in 0..10_000 {
        let mut rng = rand::thread_rng();
        let order = rng.gen_range(1..200);
        let input = buttap(order as _);
        let wo = rng.gen_range(0.0..1.0);
        let bw = rng.gen_range(0.0..1.0);
        let input2 = input.clone();
        let rust = lp2bp_zpk(input.clone(), wo, bw);

        let z_formatted = format!("{}", input2.z).replace("i", "j");
        let p_formatted = format!("{}", input2.p).replace("i", "j");

        let py_code = format!(
            "signal.lp2bp_zpk({z_formatted},{p_formatted},{}, wo={wo}, bw={bw})",
            input.k
        );

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);

        assert!(check_zpk_filter(rust, python))
    }
}

#[test]
fn with_py_fuzzy_lp2bs_zpk() {
    for _ in 0..10_000 {
        let mut rng = rand::thread_rng();
        let order = rand::thread_rng().gen_range(1..200);
        let input = buttap(order as _);
        let wo = rng.gen_range(0.0..1.0);
        let bw = rng.gen_range(0.0..1.0);

        let z_formatted = format!("{}", input.z).replace("i", "j");
        let p_formatted = format!("{}", input.p).replace("i", "j");
        let k = input.k;

        let py_code = format!(
            "signal.lp2bs_zpk({z_formatted},{p_formatted},{}, wo={wo}, bw={bw})",
            k
        );

        let rust = lp2bs_zpk(input.clone(), wo, bw);

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);

        assert!(check_zpk_filter(rust, python))
    }
}
