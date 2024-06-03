use nalgebra::RealField;
use num::traits::FloatConst;

use crate::odr::polynomial::PolynomialCoef;
use crate::optimize::*;

use crate::optimize::root_scalar::*;
use crate::optimize::util::Espilon;

pub trait BracketSolver<F, R, M>: IterativeSolver<R, R, R, R, M>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    const DEFAULT_BRACKER: (R, R);
    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self
    where
        Self: Sized;
}

#[derive(Debug, Clone)]
pub enum BracketMethod {
    Bisect,
    RegularFalsi,
    /// Ridders, C. (1979).
    /// "A new algorithm for computing a single root of a real continuous function".
    /// IEEE Transactions on Circuits and Systems. 26: 979–980. doi:10.1109/TCS.1979.1084580
    Ridder,
    Brent,
    InverseCubic,
}

pub enum BracketMethodSolver<F, R> {
    Bisect(BisectSolver<F, R>),
    RegularFalsi(RegularFalsiSolver<F, R>),
    /// Ridders, C. (1979).
    /// "A new algorithm for computing a single root of a real continuous function".
    /// IEEE Transactions on Circuits and Systems. 26: 979–980. doi:10.1109/TCS.1979.1084580
    Ridder(Ridder<F, R>),
    Brent(Brent<F, R>),
    InverseCubic(InverseCubic<F, R>),
}

impl<F, R, M> IterativeSolver<R, R, R, R, M> for BracketMethodSolver<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        match self {
            Self::Bisect(solver) => solver.new_solution(),
            Self::RegularFalsi(solver) => solver.new_solution(),
            Self::Ridder(solver) => solver.new_solution(),
            Self::Brent(solver) => solver.new_solution(),
            Self::InverseCubic(solver) => solver.new_solution(),
        }
    }
}

impl BracketMethod {
    fn valid_bracket<R, M>(bracket_x: (R, R), bracket_f: (R, R)) -> Option<String>
    where
        R: IntoMetric<M> + Float,
        M: Metric,
    {
        let (a, b) = bracket_x;
        let (fa, fb) = bracket_f;

        if a.is_infinite() || b.is_infinite() {
            return Some("Bracket must be finite".to_string());
        }
        if fa.is_infinite() || fb.is_infinite() {
            return Some("Bracket must be evaluated to be finite".to_string());
        }
        if (fa * fb).is_sign_positive() {
            return Some("Bracket must be evaluated to be different sign".to_string());
        }
        None
    }
    fn get_solver<F: FnMut(R) -> R, R, M>(
        &self,
        fun: F,
        bracket_x: (R, R),
        bracket_f: (R, R),
    ) -> BracketMethodSolver<F, R>
    where
        R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
        M: Metric,
    {
        match self {
            Self::Bisect => {
                BracketMethodSolver::Bisect(BisectSolver::new(fun, bracket_x, bracket_f))
            }
            Self::RegularFalsi => BracketMethodSolver::RegularFalsi(RegularFalsiSolver::new(
                fun, bracket_x, bracket_f,
            )),
            Self::Ridder => BracketMethodSolver::Ridder(Ridder::new(fun, bracket_x, bracket_f)),
            Self::Brent => BracketMethodSolver::Brent(Brent::new(fun, bracket_x, bracket_f)),
            Self::InverseCubic => {
                BracketMethodSolver::InverseCubic(InverseCubic::new(fun, bracket_x, bracket_f))
            }
            _ => todo!(),
        }
    }
}

pub fn solve_from_bracket<F, R, M>(
    fun: F,
    bracket_method: &BracketMethod,
    bracket: (R, R),
    criteria: Option<OptimizeCriteria<R, R, M>>,
) -> OptimizeResult<R, R, R, R, M>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: Fn(R) -> R,
{
    let evaluator = RootScalarEvaluator::new(criteria);
    let evaluator = Rc::new(RefCell::new(evaluator));

    let fun = {
        let evaluator = evaluator.clone();
        move |x| {
            evaluator.borrow_mut().res.fev();
            fun(x)
        }
    };
    let fun = Box::new(fun);

    let bracket_x = bracket;
    let bracket_f = (fun(bracket_x.0), fun(bracket_x.1));

    if let Some(msg) = BracketMethod::valid_bracket(bracket_x, bracket_f) {
        evaluator.borrow_mut().res.set_failure(msg);
    }

    let solver = bracket_method.get_solver(fun, bracket_x, bracket_f);

    iterative_optimize(solver, evaluator)
}

pub struct BisectSolver<F, R> {
    fun: F,
    pos: R,
    neg: R,
}

impl<F, R, M> BracketSolver<F, R, M> for BisectSolver<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    const DEFAULT_BRACKER: (R, R) = todo!();

    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self {
        let (pos, neg) = {
            if bracket_f.0.is_sign_positive() {
                bracket_x
            } else {
                (bracket_x.1, bracket_x.0)
            }
        };
        Self { fun, pos, neg }
    }
}
impl<F, R, M> IterativeSolver<R, R, R, R, M> for BisectSolver<F, R>
where
    R: IntoMetric<M> + Float,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        let x = (self.pos + self.neg) / R::from(2).unwrap();
        let f = (self.fun)(x);
        if f.is_sign_positive() {
            self.pos = x;
        } else {
            self.neg = x;
        }
        (x, f, None, None)
    }
}

pub fn regular_falsi<R>(pos_x: R, neg_x: R, pos_f: R, neg_f: R) -> R
where
    R: Float,
{
    (pos_f * neg_x - neg_f * pos_x) / (pos_f - neg_f + R::epsilon())
}

pub struct RegularFalsiSolver<F, R> {
    fun: F,
    pos_x: R,
    neg_x: R,
    pos_f: R,
    neg_f: R,
}

impl<F, R, M> BracketSolver<F, R, M> for RegularFalsiSolver<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    const DEFAULT_BRACKER: (R, R) = todo!();
    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self {
        let (pos_x, neg_x, pos_f, neg_f) = {
            if bracket_f.0.is_sign_positive() {
                (bracket_x.0, bracket_x.1, bracket_f.0, bracket_f.1)
            } else {
                (bracket_x.1, bracket_x.0, bracket_f.1, bracket_f.0)
            }
        };
        Self {
            fun,
            pos_x,
            neg_x,
            pos_f,
            neg_f,
        }
    }
}
impl<F, R, M> IterativeSolver<R, R, R, R, M> for RegularFalsiSolver<F, R>
where
    R: IntoMetric<M> + Float,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        let x = regular_falsi(self.pos_x, self.neg_x, self.pos_f, self.neg_f);

        let f = (self.fun)(x);

        /*
        println!("{}:{}, {}:{}, {}:{}",
        self.pos_x,self.pos_f,
        self.neg_x,self.neg_f,
        x,f
        );
        */

        if f.is_sign_positive() {
            self.pos_x = x;
            self.pos_f = f;
        } else {
            self.neg_x = x;
            self.neg_f = f;
        }

        (x, f, None, None)
    }
}

pub struct Ridder<F, R> {
    fun: F,
    x0: R,
    x2: R,
    f0: R,
    f2: R,
}

impl<F, R, M> BracketSolver<F, R, M> for Ridder<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    const DEFAULT_BRACKER: (R, R) = todo!();
    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self {
        Self {
            fun,
            x0: bracket_x.0,
            x2: bracket_x.1,
            f0: bracket_f.0,
            f2: bracket_f.1,
        }
    }
}
impl<F, R, M> IterativeSolver<R, R, R, R, M> for Ridder<F, R>
where
    R: IntoMetric<M> + Float,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        let x1 = (self.x0 + self.x2) / R::from(2.0).unwrap();
        let d = self.x0 - x1;

        let y1 = (self.fun)(x1);

        let x3 = x1
            + d * (y1 / self.f2)
                / (((y1 / self.f0).powi(2) - self.f2 / self.f0).sqrt() + R::epsilon());
        let y3 = (self.fun)(x3);

        if (y1 / y3).is_sign_negative() {
            self.x0 = x1;
            self.f0 = y1;
            self.x2 = x3;
            self.f2 = y3;
        } else if (self.f2 / y3).is_sign_negative() {
            self.x0 = x3;
            self.f0 = y3;
        } else {
            self.x2 = x3;
            self.f2 = y3;
        }

        (x3, y3, None, None)
    }
}

/// <https://mathsfromnothing.au/brents-method/?i=1>
pub struct Brent<F, R> {
    fun: F,
    a: R,
    b: R,
    c: R,
    d: Option<R>,
    fa: R,
    fb: R,
    fc: R,
    used_bisect: bool,
}

impl<F, R, M> BracketSolver<F, R, M> for Brent<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    const DEFAULT_BRACKER: (R, R) = todo!();
    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self {
        let (mut a, mut b) = bracket_x;
        let (mut fa, mut fb) = bracket_f;

        if Float::abs(fa) < Float::abs(fb) {
            (a, b) = (b, a);
            (fa, fb) = (fb, fa);
        }

        let (c, fc) = (a, fa);

        Self {
            fun,
            a,
            b,
            c,
            d: None,
            fa,
            fb,
            fc,
            used_bisect: false,
        }
    }
}
impl<F, R, M> IterativeSolver<R, R, R, R, M> for Brent<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        //println!("----");
        //println!("a: {:>8}, fa: {:>8}",self.a, self.fa);
        //println!("b: {:>8}, fb: {:>8}",self.b, self.fb);
        //println!("c: {:>8}, fb: {:>8}",self.c, self.fc);
        //println!("d: {:?}",self.d);
        //println!("u: {}", self.used_bisect);

        let mut s = if self.fa != self.fb && self.fa != self.fc && self.fb != self.fc {
            let y = [self.a, self.b, self.c];
            let x = [self.fa, self.fb, self.fc];
            least_square::poly_fit(&x, &y, 2).unwrap().eval(R::zero())
        } else {
            regular_falsi(self.a, self.b, self.fa, self.fb)
        };

        let condition1 = ((s - (R::from(3.0).unwrap() * self.a + self.b) / R::from(4.0).unwrap())
            * (s - self.b))
            .is_sign_positive();

        let condition2 = self.used_bisect
            && Float::abs(s - self.b) >= Float::abs(self.b - self.c) / R::from(2.0).unwrap();

        let condition3 = if let Some(d) = self.d {
            !self.used_bisect
                && Float::abs(s - self.b) >= Float::abs(self.c - d) / R::from(2.0).unwrap()
        } else {
            false
        };

        if condition1 || condition2 || condition3 {
            s = (self.a + self.b) / R::from(2.0).unwrap();
        }

        //println!("s : {}", s);

        let fs = (self.fun)(s);

        self.d = Some(self.c);
        self.c = self.b;
        self.fc = self.fb;

        if (self.fa * fs).is_sign_negative() {
            self.b = s;
            self.fb = fs;
        } else {
            self.a = s;
            self.fa = fs;
        }

        if Float::abs(self.fa) < Float::abs(self.fb) {
            (self.a, self.b, self.fa, self.fb) = (self.b, self.a, self.fb, self.fa);
        }

        (s, fs, None, None)
    }
}

pub struct InverseCubic<F, R> {
    fun: F,
    a: (R, R),
    b: Option<(R, R)>,
    c: Option<(R, R)>,
    d: (R, R),
}

impl<F, R, M> BracketSolver<F, R, M> for InverseCubic<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    const DEFAULT_BRACKER: (R, R) = todo!();
    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self {
        Self {
            fun,
            a: (bracket_x.0, bracket_f.0),
            b: None,
            c: None,
            d: (bracket_x.1, bracket_f.1),
        }
    }
}
impl<F, R, M> IterativeSolver<R, R, R, R, M> for InverseCubic<F, R>
where
    R: IntoMetric<M> + Float + Espilon + PolynomialCoef + RealField,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        let (a, fa) = self.a;
        let (d, fd) = self.d;

        match (self.b, self.c) {
            // For the first guess we only have the bracket,
            (None, None) => {
                // we just use regular falsi to find a new guess
                let mut g = regular_falsi(a, d, fa, fd);

                // evaluate fun at the new guess
                let fg = (self.fun)(g);

                // check the if the sign of fg is the same as fa or fd
                if (fg * fa).is_positive() {
                    // if it is the same as fa, we have b and fb, the inner guess on the a's side
                    self.b = Some((g, fg));
                } else {
                    // if it is the same as fd, we have c and fc, the inner guess on the d's side
                    self.c = Some((g, fg));
                }

                (g, fg, None, None)
            }

            // when there are no inner guess on d's side
            (Some((b, fb)), None) => {
                // we have three point, and can perform inverse quadratic interpolation
                let y = [a, b, d];
                let x = [fa, fb, fd];
                let g = least_square::poly_fit(&x, &y, 2).unwrap().eval(R::zero());

                // if the new guess does not fall in between a's inner guess and d
                if ((g - b) * (g - d)).is_positive() {
                    // we make a new guess using bisection
                    let m = (b + d) / R::from_f64(2.0).unwrap();
                    let fm = (self.fun)(m);

                    // determint if fm is on d's side
                    if (fm * fd).is_positive() {
                        // if it is, we get inner guess on d's side
                        self.c = Some((m, fm));
                    } else {
                        // if it isn't, we enclose the guess on a's side
                        self.a = (b, fb);
                        self.b = Some((m, fm));
                    }
                    // and return the new guess
                    return (m, fm, None, None);
                }

                // if the new guess does fall in between a's inner guess and d
                // we evaluat fun at the new guess
                let mut fg = (self.fun)(g);

                // if fg is on d's side
                if (fg * fd).is_positive() {
                    // we get out second guess on d's side
                    self.c = Some((g, fg));
                    return (g, fg, None, None);
                }
                // if it isn't on d's side, we use bisection to make one more guess
                let m = (b + d) / R::from_f64(2.0).unwrap();
                let fm = (self.fun)(m);

                // if fm is on d's side
                if (fm * fd).is_positive() {
                    // we get inner guess on d's side
                    self.c = Some((m, fm));
                    return (m, fm, None, None);
                }

                // if both guess are not at d's side
                // we check which one is better at enclosing
                if ((g - b) * (g - m)).is_negative() {
                    self.a = (g, fg);
                    self.b = Some((m, fm));
                    return (m, fm, None, None);
                }
                self.a = (m, fm);
                self.b = Some((g, fg));
                (g, fg, None, None)
            }

            (None, Some((c, fc))) => {
                let y = [a, c, d];
                let x = [fa, fc, fd];
                let g = least_square::poly_fit(&x, &y, 2).unwrap().eval(R::zero());

                if ((g - c) * (g - a)).is_positive() {
                    let m = (a + c) / R::from_f64(2.0).unwrap();
                    let fm = (self.fun)(m);

                    if (fm * fd).is_negative() {
                        self.d = (c, fc);
                    }
                    self.b = Some((m, fm));

                    return (m, fm, None, None);
                }

                let mut fg = (self.fun)(g);

                if (fg * fa).is_positive() {
                    self.b = Some((g, fg));
                    return (g, fg, None, None);
                }

                let m = (a + c) / R::from_f64(2.0).unwrap();
                let fm = (self.fun)(m);

                if (fm * fd).is_positive() {
                    self.c = Some((g, fm));
                    return (g, fm, None, None);
                }

                if ((g - c) * (g - m)).is_negative() {
                    self.d = (g, fg);
                    self.c = Some((m, fm));
                    return (m, fm, None, None);
                }

                self.d = (m, fm);
                self.c = Some((g, fg));

                (g, fg, None, None)
            }
            // when we have both inner and outter guess on both side
            (Some((b, fb)), Some((c, fc))) => {
                // we can use inverse cubic interpolation to generate a new guess
                let y = [a, b, c, d];
                let x = [fa, fb, fc, fd];
                let g = least_square::poly_fit(&x, &y, 3).unwrap().eval(R::zero());

                // if the new guess does not fall between inner guess of either side
                let g = if ((g - b) * (g - c)).is_positive() {
                    // change the new guess to using bisection
                    (b + c) / R::from_f64(2.0).unwrap()
                } else {
                    g
                };

                let fg = (self.fun)(g);

                if (fg * fa).is_positive() {
                    self.a = (b, fb);
                    self.b = Some((g, fg));
                } else {
                    self.d = (c, fc);
                    self.c = Some((g, fg));
                }

                (g, fg, None, None)
            }
        }
    }
}
