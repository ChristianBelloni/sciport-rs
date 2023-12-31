use crate::common::with_scipy;
use lazy_static::lazy_static;
use rand::Rng;
use sciport_rs::signal::windows::*;

lazy_static! {
    static ref TEST_LEN: u64 = std::option_env!("TEST_LEN")
        .map(str::parse)
        .and_then(Result::ok)
        .unwrap_or(1_000);
    static ref TEST_ITER: usize = std::option_env!("TEST_ITER")
        .map(str::parse)
        .and_then(Result::ok)
        .unwrap_or(1000);
}
fn len(l: u64) -> u64 {
    rand::thread_rng().gen_range(1..l)
}

fn py_bool(b: bool) -> &'static str {
    if b {
        "True"
    } else {
        "False"
    }
}

#[test]
pub fn test_boxcar() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = boxcar(len, sym).to_vec();
        let py_script = format!("signal.windows.boxcar({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();

        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_triang() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = triang(len, sym).to_vec();
        let py_script = format!("signal.windows.triang({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();

        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_blackman() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = blackman(len, sym).to_vec();
        let py_script = format!("signal.windows.blackman({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_hamming() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = hamming(len, sym).to_vec();
        let py_script = format!("signal.windows.hamming({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_hann() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = hann(len, sym).to_vec();
        let py_script = format!("signal.windows.hann({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_bartlett() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = bartlett(len, sym).to_vec();
        let py_script = format!("signal.windows.bartlett({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_flattop() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = flattop(len, sym).to_vec();
        let py_script = format!("signal.windows.flattop({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_parzen() {
    for _ in 0..*TEST_ITER {
        let len = 10;
        let sym = rand::random();
        let rust_res = parzen(len, sym).to_vec();
        let py_script = format!("signal.windows.parzen({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_bohman() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = bohman(len, sym).to_vec();
        let py_script = format!("signal.windows.bohman({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_blackmanharris() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = blackmanharris(len, sym).to_vec();
        let py_script = format!("signal.windows.blackmanharris({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_nuttall() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = nuttall(len, sym).to_vec();
        let py_script = format!("signal.windows.nuttall({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_barthann() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = barthann(len, sym).to_vec();
        let py_script = format!("signal.windows.barthann({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_cosine() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = cosine(len, sym).to_vec();
        let py_script = format!("signal.windows.cosine({len}, {})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_exponential() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = exponential(len, None, None, sym).to_vec();
        let py_script = format!("signal.windows.exponential({len}, sym={})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_tukey() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let alpha = rand::thread_rng().gen_range(0.0..1.0);
        let rust_res = tukey(len, alpha, sym).to_vec();
        let py_script = format!(
            "signal.windows.tukey({len}, alpha={alpha}, sym={})",
            py_bool(sym)
        );
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_taylor() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = taylor(len, None, None, None, sym).to_vec();
        let py_script = format!(
            "signal.windows.taylor({len},norm=False,sym={})",
            py_bool(sym)
        );
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();

        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_lanczos() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let rust_res = lanczos(len, sym).to_vec();
        let py_script = format!("signal.windows.lanczos({len}, sym={})", py_bool(sym));
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice());
    }
}

#[test]
pub fn test_kaiser() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let sym = rand::random();
        let beta = rand::thread_rng().gen_range(2.0..30.0);
        let rust_res = kaiser(len, beta, sym).to_vec();
        let py_script = format!(
            "signal.windows.kaiser({len}, beta={beta}, sym={})",
            py_bool(sym)
        );
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(
            rust_res.as_slice(),
            py_res.as_slice(),
            epsilon = 0.00000000000001
        );
    }
}

#[test]
pub fn test_kaiser_bessel_derived() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let len = if len % 2 == 1 { len + 1 } else { len };
        let sym = true;
        let beta = rand::thread_rng().gen_range(2.0..30.0);
        let rust_res = kaiser_bessel_derived(len, beta, sym).to_vec();
        let py_script = format!(
            "signal.windows.kaiser_bessel_derived({len}, beta={beta}, sym={})",
            py_bool(sym)
        );
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(
            rust_res.as_slice(),
            py_res.as_slice(),
            epsilon = 0.00000000000001
        );
    }
}

#[test]
pub fn test_gaussian() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let len = if len % 2 == 1 { len + 1 } else { len };
        let sym = true;
        let std_dev = rand::thread_rng().gen_range(2.0..30.0);
        let rust_res = gaussian(len, std_dev, sym).to_vec();
        let py_script = format!(
            "signal.windows.gaussian({len}, {std_dev}, sym={})",
            py_bool(sym)
        );
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice(),);
    }
}

#[test]
pub fn test_general_gaussian() {
    for _ in 0..*TEST_ITER {
        let len = len(*TEST_LEN);
        let len = if len % 2 == 1 { len + 1 } else { len };
        let sym = true;
        let p = rand::thread_rng().gen_range(0.0..1.0);
        let sig = rand::thread_rng().gen_range(2.0..30.0);
        let rust_res = general_gaussian(len, p, sig, sym).to_vec();
        let py_script = format!(
            "signal.windows.general_gaussian({len}, {p}, {sig}, sym={})",
            py_bool(sym)
        );
        let py_res: Vec<f64> = with_scipy(&py_script).unwrap();
        approx::assert_relative_eq!(rust_res.as_slice(), py_res.as_slice(),);
    }
}
