use crate::signal::{BandFilter, FirwinWindow};
use crate::special::sinc;
use ndarray::{Array, Array1, IxDyn};

pub fn firwin(
    numtaps: i64,
    cutoff: BandFilter,
    width: Option<f64>,
    window: FirwinWindow,
    scale: bool,
    fs: Option<f64>,
) -> Array<f64, IxDyn> {
    let nyq = 0.5 * fs.unwrap_or(2.0);
    let cutoff = cutoff / nyq;

    let pass_zero = match &cutoff {
        BandFilter::Bandpass { low: _, high: _ } | BandFilter::Highpass(_) => false,
        BandFilter::Bandstop { low: _, high: _ } | BandFilter::Lowpass(_) => true,
    };

    let pass_nyquist = ((cutoff.size() & 1) == 1) ^ pass_zero;

    if pass_nyquist && numtaps % 2 == 0 {
        panic!("A filter with an even number of coefficients must have zero response at the Nyquist frequency.");
    }

    let mut cutoff = cutoff.to_vec();
    if pass_zero {
        cutoff.insert(0, 0.0);
    }
    if pass_nyquist {
        cutoff.push(1.0);
    }

    let size = cutoff.len();
    if size % 2 != 0 {
        panic!("");
    }

    let cutoff = Array1::from_vec(cutoff);
    let bands = cutoff.into_shape((size / 2, 2)).unwrap();

    let alpha = 0.5 * ((numtaps as f64) - 1.0);

    let m = (0..numtaps)
        .map(|a| (a as f64) - alpha)
        .collect::<Array1<f64>>();

    let mut h = Array1::from_vec(vec![0.0; numtaps as usize]);

    for row in bands.rows() {
        let left = row[0];
        let right = row[1];

        h += &(right * sinc(right * m.clone()));
        h -= &(left * sinc(left * m.clone()));
    }

    // TODO!
    // implement get_window here and call it!
    //
    // h *= get_window(window, numtaps, fftbins = False);

    if scale {
        let first = bands.rows().into_iter().next().unwrap();
        let (left, right) = (first[0], first[1]);

        let scale_frequency = if left == 0.0 {
            0.0
        } else if right == 1.0 {
            1.0
        } else {
            0.5 * (left + right)
        };

        let c = (m * std::f64::consts::PI * scale_frequency).mapv(|a| a.cos());
        let s = (c * &h).sum();
        h /= s;
    }

    h.into_dyn()
}