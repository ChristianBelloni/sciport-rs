use std::cell::RefCell;
use std::rc::Rc;

use num::{complex::ComplexFloat, Float};

use crate::optimize::*;

pub mod bracket;
pub mod fixed_point;
pub mod halley;
pub mod newton;
pub mod polynomial;
pub mod secant;

pub use bracket::solve_from_bracket;
pub use fixed_point::fixed_point_method;
pub use halley::{halley_method, halley_method_approx};
pub use newton::{newton_method, newton_method_approx};
pub use polynomial::polynomial_roots;
pub use secant::secant_method;

struct RootScalarEvaluator<C, M>
where
    C: IntoMetric<M> + ComplexFloat,
    M: Metric,
{
    res: OptimizeResult<C, C, C, C, M>,
}

impl<C, M> Evaluator<C, C, C, C, M> for RootScalarEvaluator<C, M>
where
    C: IntoMetric<M> + ComplexFloat,
    M: Metric,
{
    fn new(criteria: Option<OptimizeCriteria<C, C, M>>) -> Self
    where
        Self: Sized,
    {
        let res = OptimizeResult::default()
            .set_criteria(criteria.unwrap_or_default())
            .set_target_f(Some(C::zero()));
        Self { res }
    }
    fn update(&mut self, new_solution: (C, C, Option<C>, Option<C>)) {
        self.res.update(new_solution);
    }
    fn eval(&mut self) {
        if self.res.satisfy_either() {
            self.res.set_success(format!(
                "Root finding successful, {} satisfied",
                self.res.satisfied().join(", ")
            ));
        } else if self.res.overran() {
            self.res.set_failure(format!(
                "Root finding fail, max iter {:?} reached",
                self.res.criteria.maxiter
            ));
        }
    }
    fn result(&self) -> OptimizeResult<C, C, C, C, M> {
        self.res.clone()
    }
    fn flag(&self) -> &OptimizeResultFlag {
        &self.res.flag
    }
}
