mod common;

use crate::common::check_zpk_filter;
use common::with_scipy;
use ndarray::Array1;
use ndarray_rand::RandomExt;
use num::{complex::Complex64, NumCast};
use rand::{distributions::Uniform, Rng};
use sciport_rs::signal::{
    band_filter::{lp2bp_zpk, lp2bs_zpk, lp2hp_zpk, lp2lp_zpk},
    butter::*,
    cheby1::cheb1ap,
    output_type::GenericZpk,
    tools::bilinear_zpk,
};

#[test]
fn with_py_fuzzy_lp2lp_zpk() {
    for _ in 0..1_000 {
        let mut rng = rand::thread_rng();
        let order = rng.gen_range(1..200);
        let input = cheb1ap(order as _, rng.gen_range(1.0..7.00));
        let wo = rng.gen_range(0.0..1.0);
        let input2 = input.clone();
        let rust = lp2lp_zpk(input.clone(), wo);

        let z_formatted = format!("{}", input2.z).replace('i', "j");
        let p_formatted = format!("{}", input2.p).replace('i', "j");

        let py_code = format!(
            "signal.lp2lp_zpk({z_formatted},{p_formatted},{}, wo={wo})",
            input.k
        );

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };
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

        let z_formatted = format!("{}", input2.z).replace('i', "j");
        let p_formatted = format!("{}", input2.p).replace('i', "j");

        let py_code = format!(
            "signal.lp2hp_zpk({z_formatted},{p_formatted},{}, wo={wo})",
            input.k
        );

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);

        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };
        assert!(check_zpk_filter(rust, python))
    }
}

#[test]
fn with_py_fuzzy_lp2bp_zpk() {
    for _ in 0..10_000 {
        let mut rng = rand::thread_rng();
        let order = rng.gen_range(2..100);
        let input = buttap(order as _);
        let wo = rng.gen_range(0.0..1.0);
        let bw = rng.gen_range(0.0..1.0);
        let input2 = input.clone();
        let rust = lp2bp_zpk(input.clone(), wo, bw);

        let z_formatted = format!("{}", input2.z).replace('i', "j");
        let p_formatted = format!("{}", input2.p).replace('i', "j");

        let py_code = format!(
            "signal.lp2bp_zpk({z_formatted},{p_formatted},{}, wo={wo}, bw={bw})",
            input.k
        );

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };
        let res = check_zpk_filter(rust.clone(), python.clone());
        if !res {
            println!("{}", py_code);
            println!(
                "{:?} \n{:?}",
                rust.cast_with_fn(|a| <f64 as NumCast>::from(a).unwrap()),
                python
            );
        }
        assert!(res)
    }
}

#[test]
fn with_py_fuzzy_lp2bs_zpk() {
    for _ in 0..10_000 {
        let order = rand::thread_rng().gen_range(1..5);

        let input = buttap(order as _);

        let wo = 0.5;
        let bw = 10.0;

        let z_formatted = format!("{}", input.z).replace('i', "j");
        let p_formatted = format!("{}", input.p).replace('i', "j");
        let k = input.k;

        let py_code = format!(
            "signal.lp2bs_zpk({z_formatted},{p_formatted},{}, wo={wo}, bw={bw})",
            k
        )
        .replace("+-", "-");

        println!("{py_code}");

        let rust = lp2bs_zpk(input.clone(), wo, bw);

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };

        assert!(check_zpk_filter(rust, python))
    }
}

#[test]
fn with_py_fuzzy_bilinear_zpk() {
    for _ in 0..1_000 {
        let z = Array1::<f64>::random(14, Uniform::new(-1.0, 1.0)).mapv(Into::into);
        let p = Array1::<f64>::random(14, Uniform::new(-1.0, 1.0)).mapv(Into::into);
        let k = rand::thread_rng().gen_range(0.0..5.0);

        let zpk = GenericZpk { z, p, k };
        let input = zpk.clone();
        let fs = rand::thread_rng().gen_range(10.0..500.0);
        let result = bilinear_zpk(zpk, fs);

        let z_formatted = format!("{}", input.z).replace('i', "j");
        let p_formatted = format!("{}", input.p).replace('i', "j");
        let k = input.k;

        let py_code = format!("signal.bilinear_zpk({z_formatted},{p_formatted},{k}, fs={fs})",)
            .replace("+-", "-");

        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&py_code);
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };

        assert!(check_zpk_filter(result, python));
    }
}
