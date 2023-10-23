use num::complex::Complex64;
use sciport_rs::odr::polynomial::Polynomial;

fn polynomial_equal(p1: &Polynomial<f64>, p2: &Polynomial<f64>) -> bool {
    p1.degree() == p2.degree()
        && p1
            .iter()
            .zip(p2.iter())
            .fold(true, |acc, (a, b)| acc && (a - b).abs() < 1e-12)
}

#[test]
pub fn polynomial_test() {
    let p1 = Polynomial::from(vec![-3.0, 1.0]);
    let p2 = Polynomial::from_iter(vec![-3.0, 1.0].iter());

    assert!(polynomial_equal(&p1, &p2));

    let p3 = Polynomial::from(vec![2.0, 3.0]);
    let p4 = Polynomial::from(vec![1.0, 1.0]) + Polynomial::from(vec![1.0, 2.0]);
    let p5 = Polynomial::from(vec![3.0, 4.0]) - Polynomial::from(vec![1.0, 1.0]);

    assert!(polynomial_equal(&p3, &p4));
    assert!(polynomial_equal(&p3, &p5));

    let p6 = Polynomial::from(vec![1.0, 2.0, 1.0]);
    let p7 = Polynomial::from(vec![1.0, 1.0]) * Polynomial::from(vec![1.0, 1.0]);
    let p8 = Polynomial::from(vec![2.0, 4.0, 2.0]) * 0.5;
    assert!(polynomial_equal(&p6, &p7));
    assert!(polynomial_equal(&p6, &p8));

    let p9 = Polynomial::from(vec![1.0, 2.0, 3.0]);
    let p10 = Polynomial::from(vec![1.0, 1.0, 1.0, 1.0]).differentiate();
    assert!(polynomial_equal(&p9, &p10));

    let p11 = Polynomial::from_roots_k(vec![1.0, 2.0, 3.0], 1.0);
    let p12 = Polynomial::from(vec![-1.0, 1.0])
        * Polynomial::from(vec![-2.0, 1.0])
        * Polynomial::from(vec![-3.0, 1.0]);
    assert!(polynomial_equal(&p11, &p12));

    let sol = p11.roots();
    println!("{:?}", sol);

    let roots = vec![
        Complex64::new(1.0, 2.0),
        Complex64::new(2.0, 3.0),
        Complex64::new(-1.0, -2.0),
    ];
    let k = Complex64::new(1.0, 0.0);
    let p13 = Polynomial::from_roots_k(roots, k);
    let sol = p13.roots();
    println!("{:?}", sol);
}
