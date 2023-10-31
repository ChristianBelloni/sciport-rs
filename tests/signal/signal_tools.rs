use approx::assert_relative_eq;
use ndarray::Array1;
use numpy::Complex64;
use sciport_rs::signal::{
    band_filter::BandFilter, cheby1::Cheby1Filter, output_type::DesiredFilterOutput,
    tools::generic_approx_complex_relative_slice_eq_dbg, Firwin1Filter, GenericFIRFilterSettings,
    GenericIIRFilterSettings, IIRFilterDesign, Sampling, WindowType,
};

use crate::common::with_scipy;

#[test]
fn bad_test() {
    let filter_gen = Firwin1Filter {
        settings: GenericFIRFilterSettings {
            numtaps: 8,
            cutoff: BandFilter::Lowpass(4.0),
            width: None,
            window: WindowType::Hamming,
            scale: false,
            sampling: Sampling::Digital { fs: 200.0 },
        },
    };

    let filter = filter_gen.firwin().ba();

    println!("filter: {:?}", filter);

    let signal = Array1::linspace(-1.0, 1.0, 200);
    let result = sciport_rs::signal::Filter::lfilter(&filter, signal.mapv(Into::into), None);

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

    let epsilon = 10.0_f64.powi(-7);

    assert_eq!(filtered.len(), python.len());

    println!("signal start {}", signal[0]);
    println!("signal end {:?}", signal.last());

    assert!(generic_approx_complex_relative_slice_eq_dbg(
        filtered.to_vec().as_slice(),
        python.as_slice(),
        epsilon,
        epsilon,
    ));
    // let filter = filter_gen.compute_filter(DesiredFilterOutput::Zpk).zpk();
    // let result =
    //     sciport_rs::signal::output_type::Filter::lfilter(&filter, signal.mapv(Into::into), None);

    // let filtered = result.filtered;
    // assert_relative_eq!(filtered.to_vec().as_slice(), python.as_slice());
}
