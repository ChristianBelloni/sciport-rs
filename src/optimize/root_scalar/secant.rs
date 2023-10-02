use crate::optimize::root_scalar::*;
use crate::optimize::util::*;


pub fn secant_method<C, M>(
    fun: Rc<dyn Fn(C) -> C>,
    x0: C,
    x1: C,
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

    let solver = SecantSolver::new(fun, x0, x1);
    let solver = Box::new(solver) as Box<dyn IterativeSolver<C, C, C, C, M>>;

    iterative_optimize(solver, evaluator)
}

pub struct SecantSolver<C>
where
    C: ComplexFloat,
{
    fun: Box<dyn FnMut(C) -> C>,
    x0: C,
    x1: C,
    f0: C,
    f1: C,
}

impl<C> SecantSolver<C>
where
    C: ComplexFloat + Espilon,
{
    fn new(mut fun: Box<dyn FnMut(C) -> C>, x0: C, x1: C) -> Self {
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

impl<C, M> IterativeSolver<C, C, C, C, M> for SecantSolver<C>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
{
    fn new_solution(&mut self) -> (C, C, Option<C>, Option<C>) {
        (self.x0, self.x1) = (
            self.x1,
            self.x1 - self.f1 * (self.x1 - self.x0) / (self.f1 - self.f0 + C::epsilon()),
        );
        (self.f0, self.f1) = ((self.fun)(self.x0), (self.fun)(self.x1));

        (self.x1, self.f1, None, None)
    }
}
