use num::Float;

use crate::optimize::*;

pub mod golden;

struct MinScalarEvaluator<R, M>
where
    R: IntoMetric<M> + Float,
    M: Metric,
{
    res: OptimizeResult<R, R, R, R, M>,
}

impl<R, M> Evaluator<R, R, R, R, M> for MinScalarEvaluator<R, M>
where
    R: IntoMetric<M> + Float,
    M: Metric,
{
    fn new(criteria: Option<OptimizeCriteria<R, R, M>>) -> Self
    where
        Self: Sized,
    {
        let res = OptimizeResult::default().set_criteria(criteria.unwrap_or_default());
        Self { res }
    }

    fn update(&mut self, new_solution: (R, R, Option<R>, Option<R>)) {
        self.res.update(new_solution);
    }

    fn eval(&mut self) {
        match (self.result().old_f, self.result().sol_f) {
            (Some(o), Some(s)) => {
                if o < s {
                    self.result_mut().set_failure(format!(
                        "Minimization fail, function is not strictly unimodal function"
                    ));
                    return;
                }
            }
            _ => {}
        }

        if self.res.old_x == self.res.sol_x{
            return;
        }

        if self.result().satisfy_either() {
            let satisfied = self.result_mut().satisfied().join(", ");
            self.result_mut()
                .set_success(format!("Minimization successful, {} satisfied", satisfied));
        } else if self.result().overran() {
            let maxiter = self.result_mut().criteria.maxiter;
            self.result_mut()
                .set_failure(format!("Minimization fail, max iter {:?} reached", maxiter));
        }
    }

    fn result(&self) -> &OptimizeResult<R, R, R, R, M> {
        &self.res
    }
    fn result_mut(&mut self) -> &mut OptimizeResult<R, R, R, R, M> {
        &mut self.res
    }
}
