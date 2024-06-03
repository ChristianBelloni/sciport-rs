use std::io::Read;

use crate::optimize::min_scalar::*;
use crate::optimize::util::*;

const GOLDEN_RATIO: f64 = 0.618_033_988_7;

pub fn golden_method<F, R, M>(
    fun: F,
    bracket: (R, R),
    criteria: Option<OptimizeCriteria<R, R, M>>,
) -> OptimizeResult<R, R, R, R, M>
where
    R: IntoMetric<M> + Float + Espilon,
    M: Metric,
    F: Fn(R) -> R,
{
    let evaluator = MinScalarEvaluator::new(criteria);
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

    let solver = GoldenSolver::new(fun, bracket_x, bracket_f);

    iterative_optimize(solver, evaluator)
}

pub struct GoldenSolver<F, R> {
    fun: F,
    a: R,
    x: Option<R>,
    y: Option<R>,
    b: R,
    fa: R,
    fx: Option<R>,
    fy: Option<R>,
    fb: R,
}

impl<F, R> GoldenSolver<F, R> {
    fn new(fun: F, bracket_x: (R, R), bracket_f: (R, R)) -> Self {
        Self {
            fun,
            a: bracket_x.0,
            x: None,
            y: None,
            b: bracket_x.1,
            fa: bracket_f.0,
            fx: None,
            fy: None,
            fb: bracket_f.1,
        }
    }
}

impl<F, R, M> IterativeSolver<R, R, R, R, M> for GoldenSolver<F, R>
where
    R: IntoMetric<M> + Float + Espilon,
    M: Metric,
    F: Fn(R) -> R,
{
    fn new_solution(&mut self) -> (R, R, Option<R>, Option<R>) {
        let d = (self.b - self.a) * R::from(GOLDEN_RATIO).unwrap();
        let x = self.x.unwrap_or_else(|| self.b - d);
        let fx = self.fx.unwrap_or_else(|| (self.fun)(x));

        let y = self.y.unwrap_or_else(|| self.a + d);
        let fy = self.fy.unwrap_or_else(|| (self.fun)(y));

        if fy < fx {
            self.a = x;
            self.fa = fx;
            self.x = Some(y);
            self.fx = Some(fy);
            self.y = None;
            self.fy = None;

            (y, fy, None, None)
        } else {
            self.b = y;
            self.fb = fy;
            self.y = Some(x);
            self.fy = Some(fx);
            self.x = None;
            self.fx = None;

            (x, fx, None, None)
        }
    }
}
