mod common;
use common::with_scipy;
use rand::Rng;
use sciport_rs::signal::{band_filter::BandFilter, firwin1::firwin, windows::WindowType};

#[test]
fn test_firwin() {
    for _ in 0..50_000 {
        let numtaps = rand::thread_rng().gen_range(0..50);
        let kind = rand::thread_rng().gen_range(0..4);

        let cutoff = match kind {
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

        if validate_firwin_input(&cutoff, numtaps) {
            continue;
        }

        let (wn, btype) = match &cutoff {
            BandFilter::Bandstop { low, high } => (format!("[{low}, {high}]"), "bandstop"),
            BandFilter::Bandpass { low, high } => (format!("[{low}, {high}]"), "bandpass"),
            BandFilter::Lowpass(data) => (format!("{data}"), "lowpass"),
            BandFilter::Highpass(data) => (format!("{data}"), "highpass"),
        };
        let rust_res = firwin(numtaps, cutoff, None, WindowType::Hamming, true, None).to_vec();

        let py_script = format!("signal.firwin({numtaps}, {wn}, pass_zero=\"{btype}\")");

        let py_res: Vec<f64> = with_scipy(&py_script);

        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice(), epsilon = 0.01);
    }
}

pub fn validate_firwin_input(cutoff: &BandFilter, numtaps: i64) -> bool {
    let pass_zero = cutoff.pass_zero();
    let pass_nyquist = cutoff.pass_nyquist(pass_zero);

    pass_nyquist && numtaps % 2 == 0
}
