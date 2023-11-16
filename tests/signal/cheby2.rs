use crate::common::check_zpk_filter;
use crate::common::with_scipy;
use num::complex::Complex64;
use rand::{thread_rng, Rng};
use sciport_rs::signal::{
    band_filter::BandFilter, cheby2::*, output_type::DesiredFilterOutput, GenericIIRFilterSettings,
    IIRFilterDesign, Sampling,
};

#[test]
fn with_py_test_cheby2() {
    for _ in 0..1000 {
        let order = rand::thread_rng().gen_range(0..200);
        let kind = rand::thread_rng().gen_range(0..4);
        let rp = rand::thread_rng().gen_range(0.0..10.0);
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
            0 => Sampling::Analog,
            1 => Sampling::Digital {
                fs: thread_rng().gen_range((3.0)..15.0),
            },
            _ => unreachable!(),
        };
        test_cheby2(order, band_filter, analog, rp);
    }
}

#[test]
fn test_cheb2ap() {
    for i in 0..1025 {
        println!("testing buttap order {i}");
        let rs = thread_rng().gen_range(0.0..15.0);
        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&format!(
            "signal.cheb2ap({i}, rs={rs})"
        ));
        let rust = cheb2ap(i, rs).expect("valid filter output");
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };
        assert!(check_zpk_filter(rust, python));
    }
}

fn test_cheby2(order: u32, band_filter: BandFilter, analog: Sampling, rs: f64) {
    let (wn, btype) = match &band_filter {
        BandFilter::Bandstop { low, high } => (format!("[{low}, {high}]"), "bandstop"),
        BandFilter::Bandpass { low, high } => (format!("[{low}, {high}]"), "bandpass"),
        BandFilter::Lowpass(data) => (format!("{data}"), "lowpass"),
        BandFilter::Highpass(data) => (format!("{data}"), "highpass"),
    };

    let (analog_s, fs) = match &analog {
        Sampling::Analog => ("True", "None".to_string()),
        Sampling::Digital { fs } => ("False", fs.to_string()),
    };
    let py_code = &format!(
        "signal.cheby2({order}, Wn={wn}, btype=\"{btype}\", output=\"zpk\", analog={analog_s}, fs={fs}, rs={rs})"
    );
    let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(py_code);

    let python = if let Some(p) = python {
        p
    } else {
        return;
    };

    let filter = Cheby2Filter {
        rs,
        settings: GenericIIRFilterSettings {
            order,
            band_filter,
            analog,
        },
    };

    let rust = filter
        .compute_filter(DesiredFilterOutput::Zpk)
        .expect("valid filter output")
        .zpk();

    let success = check_zpk_filter(rust.clone(), python.clone());
    if !success {
        println!("order {order} filter: {band_filter:#?}, analog {analog:#?}, rs: {rs}");

        println!("rust: {:?}", rust);
        println!("python: {:?}", python);
        println!("python code: {}", py_code);
    }
    assert!(success);
}
