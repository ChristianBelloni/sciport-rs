/*
    test for sciport-rs::optimize::root_scalar::*

    to run this test, run: cargo test root_scalar -- --nocapture
*/

use sciport_rs::optimize::root_scalar::*;

const BRACKET: (f64,f64) = (-15.0,12.0);
const ROOT: f64 = 4_f64;

fn have_root(x:f64)->f64{
    /* 
        f(x) = (x^2 + 1) * (x - 4) 
    */
    (x.powi(2) + 1_f64) * (x - ROOT)
}
fn have_root_2(x:f64)->f64{
    /* 
        f(x) = (x^2 + x + 1) * (x - 4) 
    */
    (x.powi(2) + x + 1_f64) * (x - ROOT)
}

fn touching_root(x:f64)->f64{
    /*
        f(x) = (x-3)^2
    */
    (x-3_f64).powi(2)
}

fn no_root(x:f64)->f64{
    /* 
        f(x) = x^2 + 1 
    */
    x.powi(2) + 1_f64
}

fn inf(_:f64)->f64{
    /*
        f(x) = infinity 
    */
    f64::INFINITY
}

#[test]
fn test_root_scalar() {
    let rootscalersolver = RootScalarSolver::default();

    let methods = vec![
        (BracketMethod::Bisect, "Bisect"),
        (BracketMethod::RegulaFalsi,"Regula Falsi"),
        (BracketMethod::Ridder,"Ridder")
    ];

    let funs:Vec<(Box<fn(f64) -> f64>, bool)> = vec![
        (Box::new(have_root),true),
        (Box::new(have_root_2),true),
        (Box::new(touching_root),false),
        (Box::new(no_root),false),
        (Box::new(inf),false),
    ];

    for (method, name) in methods{
        println!("[{}]",name);

        for (fun,sucess) in funs.iter(){
            let res = rootscalersolver.solve_with_bracket(&fun, BRACKET, method);
            
            println!("  {:?}",res);
            
            assert_eq!(res.flag.sucess(),*sucess);
            if res.flag.sucess(){
                assert!((res.root - ROOT).abs() < 1e-9);
            }
        }
    }
}