use crate::common::check_zpk_filter;
use crate::common::with_scipy;
use num::{complex::Complex64, NumCast};
use rand::{thread_rng, Rng};
use sciport_rs::signal::{
    band_filter::BandFilter, cheby1::*, output_type::DesiredFilterOutput, FilterDesign,
    GenericFilterSettings, Sampling,
};

#[test]
fn with_py_test_cheby1() {
    for _ in 0..10_000 {
        let order = rand::thread_rng().gen_range(0..100);

        let kind = rand::thread_rng().gen_range(0..4);
        let rp = rand::thread_rng().gen_range(0.0..10.0);
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
            0 => Sampling::Analog,
            1 => {
                let fs = thread_rng().gen_range((100.0)..500.0);
                band_filter = band_filter * fs / 2.0;
                Sampling::Digital { fs }
            }
            _ => unreachable!(),
        };
        test_cheby1(order, band_filter, analog, rp);
    }
}

#[test]
fn test_cheb1ap() {
    for i in 0..500 {
        println!("testing buttap order {i}");
        let i = thread_rng().gen_range(0..100);
        let rp = thread_rng().gen_range(0.0..20.0);
        let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(&format!(
            "signal.cheb1ap({i}, rp={rp})"
        ));
        let rust = cheb1ap(i, rp);
        let python = if let Some(p) = python {
            p
        } else {
            continue;
        };
        assert!(check_zpk_filter(rust, python));
    }
}

fn test_cheby1(order: u32, band_filter: BandFilter, analog: Sampling, rp: f64) {
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
        "signal.cheby1({order}, Wn={wn}, btype=\"{btype}\", output=\"zpk\", analog={analog_s}, fs={fs}, rp={rp})"
    );
    let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>, f64)>(py_code);
    let python = if let Some(p) = python {
        p
    } else {
        return;
    };

    let filter = Cheby1Filter {
        rp,
        settings: GenericFilterSettings {
            order,
            band_filter,
            analog,
        },
    };

    let rust = filter.filter(DesiredFilterOutput::Zpk).zpk();

    let success = check_zpk_filter(rust.clone(), python.clone());
    if !success {
        println!("order {order} filter: {band_filter:?}, analog {analog:?}, rp: {rp}");
        let rust = rust.cast_with_fn(|a| <f64 as NumCast>::from(a).unwrap());
        println!("rust: {:?}", rust);
        println!("python: {:?}", python);
        println!("python code: {}", py_code);
    }
    assert!(success);
}
