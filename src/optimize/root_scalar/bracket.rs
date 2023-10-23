use num::traits::FloatConst;

use crate::optimize::*;

use crate::optimize::root_scalar::*;

pub trait BracketSolver<F, R, M>: IterativeSolver<R, R, R, R, M>
where
    R: IntoMetric<M> + Float,
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
}

pub enum BracketMethodSolver<F, R> {
    Bisect(BisectSolver<F, R>),
    RegularFalsi(RegularFalsiSolver<F, R>),
    /// Ridders, C. (1979).
    /// "A new algorithm for computing a single root of a real continuous function".
    /// IEEE Transactions on Circuits and Systems. 26: 979–980. doi:10.1109/TCS.1979.1084580
    Ridder(Ridder<F, R>),
    Brent(Brent<F, R>),
}

impl<F, R, M> IterativeSolver<R, R, R, R, M> for BracketMethodSolver<F, R>
where
    R: IntoMetric<M> + Float,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        match self {
            BracketMethodSolver::Bisect(solver) => solver.new_solution(),
            BracketMethodSolver::RegularFalsi(solver) => solver.new_solution(),
            BracketMethodSolver::Ridder(solver) => solver.new_solution(),
            BracketMethodSolver::Brent(solver) => solver.new_solution(),
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
        R: IntoMetric<M> + Float,
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
    R: IntoMetric<M> + Float + FloatConst,
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
    R: IntoMetric<M> + Float,
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
    R: IntoMetric<M> + Float,
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
    R: IntoMetric<M> + Float,
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

/// https://mathsfromnothing.au/brents-method/?i=1
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
    R: IntoMetric<M> + Float,
    M: Metric,
    F: FnMut(R) -> R,
{
    const DEFAULT_BRACKER: (R, R) = todo!();
    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self {
        let (mut a, mut b) = bracket_x;
        let (fa, fb) = bracket_f;

        if fa.abs() < fb.abs() {
            (a, b) = (b, a);
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
    R: IntoMetric<M> + Float,
    M: Metric,
    F: FnMut(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        let mut s = if self.a != self.b && self.a != self.c {
            todo!("Inverse quadratic interpolation is not yet implemented");
        } else {
            regular_falsi(self.a, self.b, self.fa, self.fb)
        };

        let condition1 = ((s - (R::from(3.0).unwrap() * self.a + self.b) / R::from(4.0).unwrap())
            * (s - self.b))
            .is_sign_positive();

        let condition2 = self.used_bisect
            && (s - self.b).abs() >= (self.b - self.c).abs() / R::from(2.0).unwrap();

        let condition3 = if let Some(d) = self.d {
            !self.used_bisect && (s - self.b).abs() >= (self.c - d).abs() / R::from(2.0).unwrap()
        } else {
            false
        };

        if condition1 || condition2 || condition3 {
            s = (self.a + self.b) / R::from(2.0).unwrap();
        }

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

        if self.fa.abs() < self.fb.abs() {
            (self.a, self.b, self.fa, self.fb) = (self.b, self.a, self.fb, self.fa);
        }

        (s, fs, None, None)
    }
}
