mod common;

use crate::common::{check_ba_filter, check_zpk_filter};
use common::with_scipy;
use num::complex::Complex64;
use rand::{thread_rng, Rng};
use sciport_rs::signal::{
    band_filter::BandFilter,
    butter::*,
    output_type::{DesiredFilterOutput, FilterOutput},
    tools::zpk2ba,
    Analog, FilterDesign, GenericFilterSettings,
};
#[test]
fn with_py_test_butter() {
    for _ in 0..1000 {
        let order = rand::thread_rng().gen_range(0..50);
        let kind = rand::thread_rng().gen_range(0..4);

        let mut band_filter = match kind {
            0 => BandFilter::Lowpass(rand::thread_rng().gen_range((0.0)..1.0)),
            1 => BandFilter::Highpass(rand::thread_rng().gen_range((0.0)..1.0)),
            2 => {
                let x1: f64 = rand::thread_rng().gen_range((0.0)..1.0);
                let x2: f64 = rand::thread_rng().gen_range((0.0)..1.0);

                let low = x1.min(x2);
                let high = x1.max(x2);
                BandFilter::Bandpass { low, high }
            }
            3 => {
                let x1: f64 = rand::thread_rng().gen_range((0.0)..1.0);
                let x2: f64 = rand::thread_rng().gen_range((0.0)..1.0);

                let low = x1.min(x2);
                let high = x1.max(x2);
                BandFilter::Bandstop { low, high }
            }
            _ => unreachable!(),
        };

        let analog = match rand::thread_rng().gen_range(0..2) {
            0 => Analog::True,
            1 => {
                let fs = thread_rng().gen_range((3.0)..15.0);
                band_filter = band_filter * fs / 2.0;
                Analog::False { fs }
            }
            _ => unreachable!(),
        };
        test_butter(order, band_filter, analog);
    }
}

#[test]
fn test_buttap() {
    for i in 0..150 {
        println!("testing buttap order {i}");
        let python =
            with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&format!("signal.buttap({i})"));
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };
        let rust = buttap::<f64>(i);
        let rust_ba = zpk2ba(rust.clone());
        let python_ba = with_scipy::<(Vec<Complex64>, Vec<Complex64>)>(&format!(
            "signal.zpk2tf(*signal.buttap({i}))"
        ));

        let python_ba = if let Some(p) = python_ba {
            p
        } else {
            continue;
        };
        assert!(check_zpk_filter(rust, python));
        assert!(check_ba_filter(rust_ba, python_ba));
    }
}

pub fn test_butter(order: u32, band_filter: BandFilter, analog: Analog) {
    let (wn, btype) = match &band_filter {
        BandFilter::Bandstop { low, high } => (format!("[{low}, {high}]"), "bandstop"),
        BandFilter::Bandpass { low, high } => (format!("[{low}, {high}]"), "bandpass"),
        BandFilter::Lowpass(data) => (format!("{data}"), "lowpass"),
        BandFilter::Highpass(data) => (format!("{data}"), "highpass"),
    };

    let (analog_s, fs) = match &analog {
        Analog::True => ("True", "None".to_string()),
        Analog::False { fs } => ("False", fs.to_string()),
    };

    let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&format!(
        "signal.butter({order}, Wn={wn}, btype=\"{btype}\", output=\"zpk\", analog={analog_s}, fs={fs})"
    ));
    let python = if let Some(p) = python {
        p
    } else {
        return;
    };

    let filter = ButterFilter {
        settings: GenericFilterSettings {
            order,
            band_filter,
            analog,
        },
    };
    let rust = filter.filter(DesiredFilterOutput::Zpk).zpk();

    if rust.z.len() != python.0.len() {
        panic!()
    };
    if rust.p.len() != python.1.len() {
        panic!()
    };

    let success = check_zpk_filter(rust.clone(), python.clone());
    if !success {
        println!("order {order} filter: {band_filter:#?}, analog {analog:#?}");

        // println!("rust: {:?}", rust);
        // println!("python: {:?}", python);
    }

    assert!(success);

    let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>)>(&format!(
        "signal.butter({order}, Wn={wn}, btype=\"{btype}\", output=\"ba\", analog={analog_s}, fs={fs})"
    ));

    let rust = FilterOutput::get_output(
        rust,
        sciport_rs::signal::output_type::DesiredFilterOutput::Ba,
    )
    .ba();
    let python = if let Some(p) = python {
        p
    } else {
        return;
    };
    let success = check_ba_filter(rust.clone(), python.clone());
    if !success {
        println!("order {order} filter: {band_filter:#?}, analog {analog:#?}");

        println!("rust: {:?}", rust);
        println!("python: {:?}", python);
    }
    assert!(success);
}
