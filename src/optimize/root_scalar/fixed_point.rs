use crate::optimize::root_scalar::*;
use crate::optimize::util::*;

pub fn fixed_point_method<F: Fn(C) -> C, C, M>(
    fun: F,
    x0: C,
    criteria: Option<OptimizeCriteria<C, C, M>>,
) -> OptimizeResult<C, C, C, C, M>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
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

    let solver = FixedPointSolver::new(fun, x0);

    iterative_optimize(solver, evaluator)
}

pub struct FixedPointSolver<F, C>
where
    C: ComplexFloat,
    F: Fn(C) -> C,
{
    fun: F,
    x0: C,
    f0: C,
}

impl<F, C> FixedPointSolver<F, C>
where
    C: ComplexFloat + Espilon,
    F: Fn(C) -> C,
{
    fn new(mut fun: F, x0: C) -> Self {
        let f0 = fun(x0);
        Self { fun, x0, f0 }
    }
}

impl<F, C, M> IterativeSolver<C, C, C, C, M> for FixedPointSolver<F, C>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
    F: Fn(C) -> C,
{
    fn new_solution(&mut self) -> (C, C, Option<C>, Option<C>) {
        self.x0 = self.x0 - self.f0;
        self.f0 = (self.fun)(self.x0);
        (self.x0, self.f0, None, None)
    }
}
