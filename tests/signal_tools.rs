mod common;

use approx::assert_relative_eq;
use ndarray::Array1;
use numpy::Complex64;
use sciport_rs::signal::{
    band_filter::BandFilter, cheby1::Cheby1Filter, output_type::DesiredFilterOutput, Analog,
    FilterDesign, GenericFilterSettings,
};

use crate::common::with_scipy;

#[test]
fn bad_test() {
    let order = 10;
    let filter = Cheby1Filter {
        rp: 1.2,
        settings: GenericFilterSettings {
            order,
            band_filter: BandFilter::Lowpass(0.05),
            analog: Analog::False { fs: 200.0 },
        },
    }
    .filter(DesiredFilterOutput::Ba)
    .ba();

    println!("filter: {:?}", filter);

    let signal = Array1::linspace(-1.0, 1.0, 5);
    let result =
        sciport_rs::signal::output_type::Filter::lfilter(&filter, signal.mapv(Into::into), None);
    let python = with_scipy::<(Vec<Complex64>, Vec<Complex64>)>(&format!(
        "signal.butter({order}, 0.05, fs=200)"
    ));
    let _python = if let Some(p) = python {
        p
    } else {
        return;
    };
    let b_formatted = filter.b.to_string().replace('i', "j");
    let a_formatted = filter.a.to_string().replace('i', "j");
    let _zi: Array1<_> = vec![0.0; filter.a.len() - 1].into();

    let py_cmd = format!("signal.lfilter({b_formatted}, {a_formatted}, {signal})");
    let python = with_scipy::<Vec<Complex64>>(&py_cmd);
    let python = if let Some(p) = python {
        p
    } else {
        return;
    };

    // let py_cmd = format!("signal.lfilter({b_formatted}, {a_formatted}, {signal}, zi=[0,0,0])");
    // let (python, _): (Vec<Complex64>, Vec<Complex64>) = with_scipy(&py_cmd);

    let filtered = result.filtered;

    assert_relative_eq!(filtered.to_vec().as_slice(), python.as_slice());
}
