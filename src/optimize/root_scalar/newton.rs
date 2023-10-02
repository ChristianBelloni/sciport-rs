use crate::optimize::root_scalar::*;
use crate::optimize::util::*;

pub fn newton_method<C, M>(
    fun: Rc<dyn Fn(C) -> C>,
    dfun: Rc<dyn Fn(C) -> C>,
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

    let dfun = {
        let evaluator = evaluator.clone();
        move |x| {
            (*evaluator).borrow_mut().res.jev();
            dfun(x)
        }
    };
    let dfun = Box::new(dfun);

    let solver = NewtonSolver::new(fun, dfun, x0);
    let solver = Box::new(solver) as Box<dyn IterativeSolver<C, C, C, C, M>>;

    iterative_optimize(solver, evaluator)
}

pub struct NewtonSolver<C>
where
    C: ComplexFloat,
{
    fun: Box<dyn FnMut(C) -> C>,
    dfun: Box<dyn FnMut(C) -> C>,
    x0: C,
    f0: C,
    j0: C,
}

impl<C> NewtonSolver<C>
where
    C: ComplexFloat + Espilon,
{
    fn new(mut fun: Box<dyn FnMut(C) -> C>, mut dfun: Box<dyn FnMut(C) -> C>, x0: C) -> Self {
        let f0 = fun(x0);
        let j0 = dfun(x0);

        Self {
            fun,
            dfun,
            x0,
            f0,
            j0,
        }
    }
}

impl<C, M> IterativeSolver<C, C, C, C, M> for NewtonSolver<C>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
{
    fn new_solution(&mut self) -> (C, C, Option<C>, Option<C>) {
        self.x0 = self.x0 - (self.f0) / (self.j0 + C::epsilon());
        self.f0 = (self.fun)(self.x0);
        self.j0 = (self.dfun)(self.x0);

        (self.x0, self.f0, Some(self.j0), None)
    }
}
