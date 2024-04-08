use crate::optimize::root_scalar::*;
use crate::optimize::util::*;

pub fn newton_method_approx<F, C, M>(
    fun: F,
    x0: C,
    criteria: Option<OptimizeCriteria<C, C, M>>,
) -> OptimizeResult<C, C, C, C, M>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
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

    let fun = Rc::new(fun);

    let dfun = {
        let evaluator = evaluator.clone();
        let f = fun.clone();

        approx_derivative(move |x| f(x))
    };

    let fun = move |x| fun(x);

    let solver = NewtonSolver::new(fun, dfun, x0);

    iterative_optimize(solver, evaluator)
}

pub fn newton_method<F, FD, C, M>(
    fun: F,
    dfun: FD,
    x0: C,
    criteria: Option<OptimizeCriteria<C, C, M>>,
) -> OptimizeResult<C, C, C, C, M>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
    F: Fn(C) -> C,
    FD: Fn(C) -> C,
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

    let dfun = {
        let evaluator = evaluator.clone();
        move |x| {
            evaluator.borrow_mut().res.jev();
            dfun(x)
        }
    };

    let solver = NewtonSolver::new(fun, dfun, x0);

    iterative_optimize(solver, evaluator)
}

pub struct NewtonSolver<F, FD, C>
where
    C: ComplexFloat + Espilon,
{
    fun: F,
    dfun: FD,
    x0: C,
    f0: C,
    j0: C,
}

impl<F, FD, C> NewtonSolver<F, FD, C>
where
    C: ComplexFloat + Espilon,
    F: Fn(C) -> C,
    FD: Fn(C) -> C,
{
    fn new(mut fun: F, mut dfun: FD, x0: C) -> Self {
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

impl<F, FD, C, M> IterativeSolver<C, C, C, C, M> for NewtonSolver<F, FD, C>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
    F: Fn(C) -> C,
    FD: Fn(C) -> C,
{
    fn new_solution(&mut self) -> (C, C, Option<C>, Option<C>) {
        self.x0 = self.x0 - (self.f0) / (self.j0 + C::epsilon());
        self.f0 = (self.fun)(self.x0);
        self.j0 = (self.dfun)(self.x0);

        (self.x0, self.f0, Some(self.j0), None)
    }
}
