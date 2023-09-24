pub fn kaiser_beta(input: f64) -> f64 {
    if input > 50.0 {
        0.1102 * (input - 8.7)
    } else if input > 21.0 {
        0.5842 * (input - 21.0).powf(0.4) + 0.07886 * (input - 21.0)
    } else {
        0.0
    }
}

pub fn kaiser_atten(numtaps: i64, width: f64) -> f64 {
    2.285 * ((numtaps - 1) as f64) * std::f64::consts::PI * width + 7.95
}

pub fn kaiserord(ripple: f64, width: f64) -> (i64, f64) {
    let ripple = ripple.abs();

    if ripple < 8.0 {
        panic!("ripple attenuation too small!")
    }

    let beta = kaiser_beta(ripple);

    let numtaps = (ripple - 7.95) / 2.285 / (std::f64::consts::PI * width) + 1.0;

    let numtaps = f64::ceil(numtaps) as i64;

    (numtaps, beta)
}
