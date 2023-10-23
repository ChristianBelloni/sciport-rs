mod common;

use crate::common::check_zpk_filter;
use common::with_scipy;
use num::complex::Complex64;
use rand::{thread_rng, Rng};
use sciport_rs::signal::{
    band_filter::BandFilter, bessel::*, output_type::DesiredFilterOutput, Analog, FilterDesign,
    GenericFilterSettings,
};

#[test]
fn with_py_test_bessel() {
    for _ in 0..500 {
        let order = rand::thread_rng().gen_range(0..53);
        let kind = rand::thread_rng().gen_range(0..4);

        let band_filter = match kind {
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
            1 => Analog::False {
                fs: thread_rng().gen_range((3.0)..15.0),
            },
            _ => unreachable!(),
        };
        test_bessel(order, band_filter, analog);
    }
}

#[test]
fn test_besselap() {
    for i in 1..53 {
        println!("testing besselap order {i}");
        let python =
            with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&format!("signal.besselap({i})"));
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };
        let rust = besselap::<f64>(i, BesselNorm::Phase);
        assert!(check_zpk_filter(rust, python));
    }
}

fn test_bessel(order: u32, band_filter: BandFilter, analog: Analog) {
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
        "signal.bessel({order}, Wn={wn}, btype=\"{btype}\", output=\"zpk\", analog={analog_s}, fs={fs})"
    ));

    let python = if let Some(p) = python {
        p
    } else {
        return;
    };

    let filter = BesselFilter {
        norm: BesselNorm::Phase,
        settings: GenericFilterSettings {
            order,
            band_filter,
            analog,
        },
    };

    let rust = filter.filter(DesiredFilterOutput::Zpk).zpk();
    let success = check_zpk_filter(rust.clone(), python.clone());
    if !success {
        println!("order {order} filter: {band_filter:#?}, analog {analog:#?}");

        //println!("rust: {:#?}", rust);
        //println!("python: {:#?}", python);
    }
    assert!(success);
}
