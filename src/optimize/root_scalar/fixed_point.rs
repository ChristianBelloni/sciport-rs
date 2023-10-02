use crate::optimize::root_scalar::*;
use crate::optimize::util::*;
use crate::optimize::*;

pub fn fixed_point_method<C, M>(
    fun: Rc<dyn Fn(C) -> C>,
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
            (*evaluator).borrow_mut().res.fev();
            fun(x)
        }
    };
    let fun = Box::new(fun);

    let solver = FixedPointSolver::new(fun, x0);
    let solver = Box::new(solver) as Box<dyn IterativeSolver<C, C, C, C, M>>;

    iterative_optimize(solver, evaluator)
}

pub struct FixedPointSolver<C>
where
    C: ComplexFloat,
{
    fun: Box<dyn FnMut(C) -> C>,
    x0: C,
    f0: C,
}

impl<C> FixedPointSolver<C>
where
    C: ComplexFloat + Espilon,
{
    fn new(mut fun: Box<dyn FnMut(C) -> C>, x0: C) -> Self {
        let f0 = fun(x0);
        Self { fun, x0, f0 }
    }
}

impl<C, M> IterativeSolver<C, C, C, C, M> for FixedPointSolver<C>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
{
    fn new_solution(&mut self) -> (C, C, Option<C>, Option<C>) {
        self.x0 = self.x0 - self.f0;
        self.f0 = (self.fun)(self.x0);
        (self.x0, self.f0, None, None)
    }
}
