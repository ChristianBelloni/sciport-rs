use itertools::{EitherOrBoth, Itertools};
use num::complex::{Complex32, Complex64, ComplexFloat};
use std::fmt::{Debug, Display};
use std::ops::{Add, Div, Index, Mul, Sub};
use std::rc::Rc;

use crate::optimize::least_square;
use crate::optimize::root_scalar::polynomial::{polynomial_roots, IntoComplex};
use crate::optimize::util::Espilon;
use crate::optimize::{IntoMetric, Metric};

/// # `PolynomialCoef`
/// a common trait for polynomial coefficient, just for display purposes
pub trait PolynomialCoef: ComplexFloat + Clone + Debug + 'static {
    fn coef_to_string(&self) -> String;
}

/// # `Polynomial`
/// a polynomial struct that represented by a `Vec` of coefficient,
/// which can be `f32`,`f64`,`Complex32`, `Complex64`
///
/// where the i-th coefficient represent the i-th power's coefficient of the polynomial
#[derive(Debug, Clone)]
pub struct Polynomial<T>
where
    T: PolynomialCoef,
{
    coef: Vec<T>,
}

impl<T> Polynomial<T>
where
    T: PolynomialCoef,
{
    /// while the highest power of the polynomial is zero,
    /// pop that coefficient of the polynomial
    ///
    /// if the polynomial ended up to have no coefficient,
    /// push a zero to represent a zero constant
    pub fn saturate(mut self) -> Self {
        while let Some(&c) = self.coef.last() {
            if c == T::zero() {
                self.coef.pop();
            } else {
                break;
            }
        }
        if self.degree() == 0 {
            self.coef.push(T::zero());
        }
        self
    }
    /// return the degree of the polynomial,
    /// aka the highest power the polynomial consisted of
    pub fn degree(&self) -> usize {
        self.coef.len().clamp(1, usize::MAX) - 1
    }
    /// iterate its coefficient, from 0-th power.
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.coef.iter()
    }
    /// construct a `Polynomial<T>` for a Vec<T>``
    pub fn from_vec(coef: Vec<T>) -> Self {
        Self { coef }.saturate()
    }
    /// return new polynomial with only a zero constant
    pub fn zero() -> Self {
        Self::from_vec(vec![T::zero()])
    }
    /// return new polynomial with only a one constant
    pub fn one() -> Self {
        Self::from_vec(vec![T::one()])
    }
    /// evaluate the polynomial at `x`
    pub fn eval(&self, x: T) -> T {
        self.iter().rev().fold(T::zero(), |acc, &c| acc * x + c)
    }
    /// evaluate the polynomial at `x` for `x` in `xs`
    pub fn eval_iter(&self, xs: impl IntoIterator<Item = T>) -> Vec<T> {
        xs.into_iter().map(|x| self.eval(x)).collect()
    }
    /// return the multiply of polynomial by `x^p`
    pub fn mul_power(&self, p: usize) -> Self {
        vec![T::zero(); p].iter().chain(self.iter()).collect()
    }
    /// return differentiated polynomial
    pub fn differentiate(&self) -> Self {
        self.iter()
            .enumerate()
            .filter_map(|(i, &c)| {
                if i == 0 {
                    None
                } else {
                    Some(T::from(i as i64).unwrap() * c)
                }
            })
            .collect()
    }
    /// construct a new polynomial from roots and multiply constant `k`
    pub fn from_roots_k(roots: impl IntoIterator<Item = T>, k: T) -> Self {
        roots
            .into_iter()
            .map(|r| Polynomial::from(vec![-r, T::one()]))
            .fold(Polynomial::one(), |acc, p| acc * p)
            * k
    }
    /// return the deflated polynomial using horner's method
    ///
    /// it return the quotient polynomial and the remainder scalar
    ///
    /// > https://en.wikipedia.org/wiki/Horner%27s_method
    pub fn deflate(&self, x: T) -> Option<(Polynomial<T>, T)> {
        let result = self
            .iter()
            .rev()
            .scan(T::zero(), |carry, &coef| {
                let new_coef = coef + *carry;
                *carry = new_coef * x;
                Some(new_coef)
            })
            .collect::<Vec<_>>();
        let (remainder, quotient) = result.split_last()?;
        Some((quotient.iter().rev().collect(), remainder.to_owned()))
    }
    /// take ownership and package the polynomial into `Rc<dyn Fn(T)->T>`
    pub fn as_rc(self) -> Rc<dyn Fn(T) -> T> {
        Rc::new(move |x| self.eval(x))
    }
    /// find all root of the polynomial,
    ///
    /// where all its root will be in complex number data structure
    /// i.e. `Complex32` or `Complexf64`
    pub fn roots<C, M>(&self) -> Vec<C>
    where
        T: PolynomialCoef + Espilon + IntoMetric<M> + IntoComplex<C>,
        C: PolynomialCoef + From<T> + Espilon + IntoMetric<M> + IntoComplex<C>,
        M: Metric,
    {
        polynomial_roots(&self)
    }
    /// calculate the polynomial least square curve fit on data `x` and `y`
    /// see `sciport_rs::optimize::least_square::poly_fit`
    pub fn poly_fit<'a, Q>(
        x: impl IntoIterator<Item = &'a T>,
        y: impl IntoIterator<Item = &'a T>,
        order: usize,
    ) -> Result<Self, String>
    where
        T: Debug + Display + ComplexField<RealField = Q> + PolynomialCoef + Espilon,
        Q: ComplexFloat + Scalar + Debug + Espilon,
    {
        least_square::poly_fit(x, y, order)
    }
}

impl<T> Index<usize> for Polynomial<T>
where
    T: PolynomialCoef,
{
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.coef[index]
    }
}

impl<T> Mul<T> for Polynomial<T>
where
    T: PolynomialCoef,
{
    type Output = Self;
    fn mul(self, rhs: T) -> Self::Output {
        self.iter().map(|&c| c * rhs).collect()
    }
}
impl<T> Div<T> for Polynomial<T>
where
    T: PolynomialCoef,
{
    type Output = Self;
    fn div(self, rhs: T) -> Self::Output {
        self * rhs.recip()
    }
}

impl<T> Add for Polynomial<T>
where
    T: PolynomialCoef,
{
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        self.iter()
            .zip_longest(rhs.iter())
            .map(|pair| match pair {
                EitherOrBoth::Both(&a, &b) => a + b,
                EitherOrBoth::Left(&a) => a,
                EitherOrBoth::Right(&b) => b,
            })
            .collect()
    }
}

impl<T> Mul for Polynomial<T>
where
    T: PolynomialCoef,
{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        self.iter()
            .enumerate()
            .map(|(i, &c)| rhs.mul_power(i) * c)
            .fold(Polynomial::zero(), |acc, p| acc + p)
    }
}

impl<T> Sub for Polynomial<T>
where
    T: PolynomialCoef,
{
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        self + rhs * (T::zero() - T::one())
    }
}

impl<T> Display for Polynomial<T>
where
    T: PolynomialCoef,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.iter()
                .enumerate()
                .map(|(i, &c)| {
                    format!(
                        "{}{}",
                        c.coef_to_string(),
                        if i == 0 {
                            format!("")
                        } else {
                            format!(" * x^{:<3}", i)
                        }
                    )
                })
                .collect::<Vec<_>>()
                .join(" + ")
        )
    }
}

impl<T> FromIterator<T> for Polynomial<T>
where
    T: PolynomialCoef,
{
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        Self {
            coef: Vec::from_iter(iter),
        }
        .saturate()
    }
}

impl<'a, T> FromIterator<&'a T> for Polynomial<T>
where
    T: PolynomialCoef,
{
    fn from_iter<I: IntoIterator<Item = &'a T>>(iter: I) -> Self {
        iter.into_iter().cloned().collect()
    }
}

impl<T> From<&[T]> for Polynomial<T>
where
    T: PolynomialCoef,
{
    fn from(value: &[T]) -> Self {
        value.into_iter().collect()
    }
}

impl<T> From<Vec<T>> for Polynomial<T>
where
    T: PolynomialCoef,
{
    fn from(value: Vec<T>) -> Self {
        Self::from_vec(value)
    }
}

impl<T> IntoIterator for Polynomial<T>
where
    T: PolynomialCoef,
{
    type Item = T;
    type IntoIter = std::vec::IntoIter<T>;
    fn into_iter(self) -> Self::IntoIter {
        self.coef.into_iter()
    }
}

impl PolynomialCoef for f64 {
    fn coef_to_string(&self) -> String {
        format!("{:9.2}", self)
    }
}
impl PolynomialCoef for f32 {
    fn coef_to_string(&self) -> String {
        format!("{:9.2}", self)
    }
}
impl PolynomialCoef for Complex32 {
    fn coef_to_string(&self) -> String {
        format!("({:9.2} + {:9.2}i)", self.re, self.im)
    }
}
impl PolynomialCoef for Complex64 {
    fn coef_to_string(&self) -> String {
        format!("({:9.2} + {:9.2}i)", self.re, self.im)
    }
}
