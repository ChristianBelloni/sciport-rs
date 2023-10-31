use crate::optimize::root_scalar::*;
use num::traits::Float;

pub fn secant_method<F, C, M>(
    fun: F,
    x0: C,
    x1: C,
    criteria: Option<OptimizeCriteria<C, C, M>>,
) -> OptimizeResult<C, C, C, C, M>
where
    C: IntoMetric<M> + ComplexFloat,
    M: Metric,
    F: Fn(C) -> C,
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

    let solver = SecantSolver::new(fun, x0, x1);

    iterative_optimize(solver, evaluator)
}

pub struct SecantSolver<F, C>
where
    C: ComplexFloat,
{
    fun: F,
    x0: C,
    x1: C,
    f0: C,
    f1: C,
}

impl<F, C> SecantSolver<F, C>
where
    C: ComplexFloat,
    F: Fn(C) -> C,
{
    fn new(mut fun: F, x0: C, x1: C) -> Self {
        let (f0, f1) = (fun(x0), fun(x1));
        Self {
            fun,
            x0,
            x1,
            f0,
            f1,
        }
    }
}

impl<F, C, M> IterativeSolver<C, C, C, C, M> for SecantSolver<F, C>
where
    C: IntoMetric<M> + ComplexFloat,
    M: Metric,
    F: Fn(C) -> C,
{
    fn new_solution(&mut self) -> (C, C, Option<C>, Option<C>) {
        (self.x0, self.x1) = (
            self.x1,
            self.x1
                - self.f1 * (self.x1 - self.x0)
                    / (self.f1 - self.f0 + C::from(C::Real::epsilon()).unwrap()),
        );
        (self.f0, self.f1) = ((self.fun)(self.x0), (self.fun)(self.x1));

        (self.x1, self.f1, None, None)
    }
}
