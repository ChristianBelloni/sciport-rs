use sciport_rs::optimize::min_scalar;
use sciport_rs::optimize::*;

#[test]
fn min_scalar() {
    let fun = |x: f64| x.powi(2);

    let bracket = (-20f64, 20f64);
    let criteria = Some(
        OptimizeCriteria::empty()
            .set_xatol(Some(1e-9f64))
            .set_maxiter(Some(5000)),
    );

    let res = min_scalar::golden::golden_method(fun, bracket, criteria);

    println!("{}", res);
}
