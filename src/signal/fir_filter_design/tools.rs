use num::{Float, NumCast};

pub fn kaiser_beta<T: Float>(input: T) -> T {
    let input = <f64 as NumCast>::from(input).unwrap();
    T::from(if input > 50.0 {
        0.1102 * (input - 8.7)
    } else if input > 21.0 {
        0.5842 * (input - 21.0).powf(0.4) + 0.07886 * (input - 21.0)
    } else {
        0.0
    })
    .unwrap()
}

pub fn kaiser_atten<T: Float>(numtaps: i64, width: T) -> T {
    let width = <f64 as NumCast>::from(width).unwrap();
    T::from(2.285 * ((numtaps - 1) as f64) * std::f64::consts::PI * width + 7.95).unwrap()
}

pub fn kaiserord<T: Float>(ripple: T, width: T) -> (i64, T) {
    let ripple = ripple.abs();

    if ripple < T::from(8.0).unwrap() {
        panic!("ripple attenuation too small!")
    }

    let beta = kaiser_beta(ripple);
    let ripple = <f64 as NumCast>::from(ripple).unwrap();
    let width = <f64 as NumCast>::from(width).unwrap();
    let numtaps = (ripple - 7.95) / 2.285 / (std::f64::consts::PI * width) + 1.0;

    let numtaps = f64::ceil(numtaps) as i64;

    (numtaps, beta)
}
