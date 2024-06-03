use std::rc::Rc;

use crate::optimize::root_scalar::*;
use crate::optimize::util::*;

pub fn halley_method_approx<F, C, M>(
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
        let f = fun.clone();

        approx_derivative(move |x| f(x))
    };

    let ddfun = {
        let f = fun.clone();

        approx_second_derivative(move |x| f(x))
    };

    let fun = move |x| fun(x);

    let solver = NewtonSolver::new(fun, dfun, ddfun, x0);

    iterative_optimize(solver, evaluator)
}

pub fn halley_method<F, FD1, FD2, C, M>(
    fun: F,
    dfun: FD1,
    ddfun: FD2,
    x0: C,
    criteria: Option<OptimizeCriteria<C, C, M>>,
) -> OptimizeResult<C, C, C, C, M>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
    F: Fn(C) -> C,
    FD1: Fn(C) -> C,
    FD2: Fn(C) -> C,
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

    let ddfun = {
        let evaluator = evaluator.clone();
        move |x| {
            evaluator.borrow_mut().res.hev();
            ddfun(x)
        }
    };

    let solver = NewtonSolver::new(fun, dfun, ddfun, x0);

    iterative_optimize(solver, evaluator)
}

pub struct NewtonSolver<F, FD1, FD2, C>
where
    C: ComplexFloat,
{
    fun: F,
    dfun: FD1,
    ddfun: FD2,
    x0: C,
    f0: C,
    j0: C,
    h0: C,
}

impl<F, FD1, FD2, C> NewtonSolver<F, FD1, FD2, C>
where
    C: ComplexFloat + Espilon,
    F: Fn(C) -> C,
    FD1: Fn(C) -> C,
    FD2: Fn(C) -> C,
{
    fn new(mut fun: F, mut dfun: FD1, mut ddfun: FD2, x0: C) -> Self {
        let f0 = fun(x0);
        let j0 = dfun(x0);
        let h0 = ddfun(x0);

        Self {
            fun,
            dfun,
            ddfun,
            x0,
            f0,
            j0,
            h0,
        }
    }
}

impl<F, FD1, FD2, C, M> IterativeSolver<C, C, C, C, M> for NewtonSolver<F, FD1, FD2, C>
where
    C: IntoMetric<M> + ComplexFloat + Espilon,
    M: Metric,
    F: Fn(C) -> C,
    FD1: Fn(C) -> C,
    FD2: Fn(C) -> C,
{
    fn new_solution(&mut self) -> (C, C, Option<C>, Option<C>) {
        self.x0 = self.x0
            - (C::from(2.0).unwrap() * self.f0 * self.j0)
                / (C::from(2.0).unwrap() * self.j0.powi(2) - self.f0 * self.h0 + C::epsilon());
        self.f0 = (self.fun)(self.x0);
        self.j0 = (self.dfun)(self.x0);
        self.h0 = (self.ddfun)(self.x0);

        (self.x0, self.f0, Some(self.j0), Some(self.h0))
    }
}
