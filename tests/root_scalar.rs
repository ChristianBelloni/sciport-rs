use num::complex::Complex64;
use std::rc::Rc;

use sciport_rs::optimize::root_scalar::{
    bracket::BracketMethod, fixed_point_method, halley_method, newton_method, secant_method,
    solve_from_bracket,
};
use sciport_rs::optimize::OptimizeCriteria;

fn print_divider(s: String) {
    println!("{}", "-".to_string().repeat(64));
    println!("| {:^60} |", s);
    println!("{}", "-".to_string().repeat(64));
}

#[test]
fn root_scalar() {
    let methods = [
        BracketMethod::Bisect,
        BracketMethod::RegularFalsi,
        BracketMethod::Ridder,
    ];
    let fun = (|x: f64| 4.0 * x.powi(3) - 12.0 * x.powi(2) + 3.0 * x - 9.0);
    let dfun = (|x: f64| 12.0 * x.powi(2) - 24.0 * x + 3.0);
    let ddfun = (|x: f64| 24.0 * x - 24.0);
    let x0 = 10.0;
    let x1 = -50.0;
    let bracket = (-10f64, 20f64);
    let criteria = Some(
        OptimizeCriteria::empty()
            .set_fltol(Some(1e-9f64))
            .set_maxiter(Some(5000)),
    );

    for bracket_method in methods {
        print_divider(format!("{:?}", bracket_method));
        let res = solve_from_bracket(fun.clone(), &bracket_method, bracket, criteria.clone());
        println!("{}", res);
    }

    print_divider("Secant".to_string());
    let res = secant_method(fun.clone(), x0, x1, criteria.clone());
    println!("{}", res);

    print_divider("Newton".to_string());
    let res = newton_method(fun.clone(), dfun.clone(), x0, criteria.clone());
    println!("{}", res);

    print_divider("Halley".to_string());
    let res = halley_method(
        fun.clone(),
        dfun.clone(),
        ddfun.clone(),
        x0,
        criteria.clone(),
    );
    println!("{}", res);

    let fun = (|z: Complex64| 2.0 * z.powi(3) - 10.0 * z.powi(2) + 3.0 * z - 15.0);
    let dfun = (|z: Complex64| 6.0 * z.powi(2) - 20.0 * z + 3.0);
    let ddfun = (|z: Complex64| 12.0 * z - 20.0);
    let x0 = Complex64::new(5.0, 10.0);
    let x1 = Complex64::new(-10.0, -20.0);
    let criteria = OptimizeCriteria::empty()
        .set_fltol(Some(1e-12f64))
        .set_maxiter(Some(1000));
    let criteria = Some(criteria);

    print_divider("Secant Complex".to_string());
    let res = secant_method(fun.clone(), x0, x1, criteria.clone());
    println!("{}", res);

    print_divider("Newton Complex".to_string());
    let res = newton_method(fun.clone(), dfun.clone(), x0, criteria.clone());
    println!("{}", res);

    print_divider("Halley Complex".to_string());
    let res = halley_method(
        fun.clone(),
        dfun.clone(),
        ddfun.clone(),
        x0,
        criteria.clone(),
    );
    println!("{}", res);

    let fun = |x: f64| x.sin();
    let x0 = 1.0;
    let criteria = Some(
        OptimizeCriteria::empty()
            .set_fltol(Some(1e-9f64))
            .set_maxiter(Some(5000)),
    );

    print_divider("Fixed Point".to_string());
    let res = fixed_point_method(fun.clone(), x0, criteria.clone());
    println!("{}", res);

    let fun = (|x: Complex64| x.sin());
    let x0 = Complex64::new(0.5, 0.3);
    let criteria = Some(
        OptimizeCriteria::empty()
            .set_fltol(Some(1e-9f64))
            .set_maxiter(Some(5000)),
    );

    print_divider("Fixed Point Complex".to_string());
    let res = fixed_point_method(fun.clone(), x0, criteria.clone());
    println!("{}", res);
}
