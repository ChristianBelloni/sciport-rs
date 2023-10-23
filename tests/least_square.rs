use numpy::Complex64;
use sciport_rs::odr::polynomial::Polynomial;
use sciport_rs::optimize::least_square;

#[test]
fn least_square_poly_fit() -> Result<(), String> {
    println!("[ Real Example ]");
    let x = (0..100).map(|v| v as f64).collect::<Vec<_>>();
    let y = x
        .iter()
        .map(|f| 10.0 - 0.25 * f + 0.11 * f * f)
        .collect::<Vec<_>>();

    let res = least_square::poly_fit(&x, &y, 1)?;
    println!("{:?}", res);
    let res = least_square::poly_fit(&x, &y, 2)?;
    println!("{:?}", res);
    let res = least_square::poly_fit(&x, &y, 3)?;
    println!("{:?}", res);

    println!("[ Complex Example ]");
    let x = (0..100)
        .map(|v| Complex64::new(v as f64, v as f64 * 2.0))
        .collect::<Vec<_>>();
    let a = Complex64::new(10.0, 5.5);
    let b = Complex64::new(-2.0, 3.3);
    let c = Complex64::new(0.22, 3.9);
    let d = Complex64::new(12.0, -9.0);
    let y = x
        .iter()
        .map(|f| d - c * f + b * f.powi(2) + a * f.powi(3))
        .collect::<Vec<_>>();
    let res = least_square::poly_fit(&x, &y, 1)?;
    println!("{:#?}", res);
    let res = least_square::poly_fit(&x, &y, 2)?;
    println!("{:#?}", res);
    let res = least_square::poly_fit(&x, &y, 3)?;
    println!("{:#?}", res);
    let res = least_square::poly_fit(&x, &y, 4)?;
    println!("{:#?}", res);

    Ok(())
}
