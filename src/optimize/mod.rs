use std::cell::RefCell;
use std::fmt::Display;
use std::rc::Rc;

pub mod criteria;
pub mod least_square;
pub mod metric;
pub mod root_scalar;
pub mod util;

pub use criteria::*;
pub use least_square::*;
pub use metric::*;

/// Iterative Optimize
pub fn iterative_optimize<S, E, X, F, J, H, M>(
    mut solver: S,
    evaluator: Rc<RefCell<E>>,
) -> OptimizeResult<X, F, J, H, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    J: IntoMetric<M>,
    H: IntoMetric<M>,
    M: Metric,
    S: IterativeSolver<X, F, J, H, M>,
    E: Evaluator<X, F, J, H, M>,
{
    while evaluator.borrow().flag().is_running() {
        let s = solver.new_solution();
        evaluator.borrow_mut().update(s);
        evaluator.borrow_mut().eval();
    }
    evaluator.borrow().result()
}

/// # Iterative Solver
/// the iterative solver that continue gaving out new solution for evaluator.
pub trait IterativeSolver<X, F, J, H, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    J: IntoMetric<M>,
    H: IntoMetric<M>,
    M: Metric,
{
    fn new_solution(&mut self) -> (X, F, Option<J>, Option<H>);
}

/// # Evaluator
/// evaluate solution given by `IterativeSolver`
pub trait Evaluator<X, F, J, H, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    J: IntoMetric<M>,
    H: IntoMetric<M>,
    M: Metric,
{
    /// create a new evaluator with a criteria
    fn new(criteria: Option<OptimizeCriteria<X, F, M>>) -> Self
    where
        Self: Sized;
    /// update self from the newsolution
    fn update(&mut self, new_solution: (X, F, Option<J>, Option<H>));
    /// evaluate, change self.flag
    fn eval(&mut self);
    /// get the result
    fn result(&self) -> OptimizeResult<X, F, J, H, M>;
    /// get the flag
    fn flag(&self) -> &OptimizeResultFlag;
}

/// # OptimizeResultFlag
/// indicate if the result is success, failure or still running
#[derive(Debug, Clone)]
pub enum OptimizeResultFlag {
    Running,
    Success(String),
    Failure(String),
}

impl OptimizeResultFlag {
    fn is_running(&self) -> bool {
        match self {
            Self::Running => true,
            _ => false,
        }
    }
    fn is_success(&self) -> bool {
        match self {
            Self::Success(_) => true,
            _ => false,
        }
    }
    fn is_failure(&self) -> bool {
        match self {
            Self::Failure(_) => true,
            _ => false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct OptimizeResult<X, F, J, H, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    J: IntoMetric<M>,
    H: IntoMetric<M>,
    M: Metric,
{
    pub criteria: OptimizeCriteria<X, F, M>,
    pub old_x: Option<X>,
    pub old_f: Option<F>,
    pub sol_x: Option<X>,
    pub sol_f: Option<F>,
    pub sol_j: Option<J>,
    pub sol_h: Option<H>,
    pub target_f: Option<F>,
    pub iter: u64,
    pub nfev: u64,
    pub njev: u64,
    pub nhev: u64,
    pub err_xa: Option<M>,
    pub err_xr: Option<M>,
    pub err_fa: Option<M>,
    pub err_fr: Option<M>,
    pub err_fl: Option<M>,
    pub flag: OptimizeResultFlag,
}

impl<X, F, J, H, M> Default for OptimizeResult<X, F, J, H, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    J: IntoMetric<M>,
    H: IntoMetric<M>,
    M: Metric,
{
    /// A default new result with all field initalizedm and default criteria
    fn default() -> Self {
        OptimizeResult {
            criteria: OptimizeCriteria::default(),
            old_x: None,
            old_f: None,
            sol_x: None,
            sol_f: None,
            sol_j: None,
            sol_h: None,
            target_f: None,
            iter: 0,
            nfev: 0,
            njev: 0,
            nhev: 0,
            err_xa: None,
            err_xr: None,
            err_fa: None,
            err_fr: None,
            err_fl: None,
            flag: OptimizeResultFlag::Running,
        }
    }
}

impl<X, F, J, H, M> Display for OptimizeResult<X, F, J, H, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    J: IntoMetric<M>,
    H: IntoMetric<M>,
    M: Metric,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // f.write_fmt(format_args!(
        //     "{}",
        //     [
        //         format!("[Optimization Result]"),
        //         format!("    flag     : {:?}", self.flag),
        //         format!("    sol_x    : {:?}", self.sol_x),
        //         format!("    sol_f    : {:?}", self.sol_f),
        //         format!("    sol_j    : {:?}", self.sol_j),
        //         format!("    sol_h    : {:?}", self.sol_h),
        //         format!("    iter     : {:?}", self.iter),
        //         format!("    target_f : {:?}", self.target_f),
        //         format!("    nfev     : {:?}", self.nfev),
        //         format!("    njev     : {:?}", self.njev),
        //         format!("    nhev     : {:?}", self.nhev),
        //         format!("    [criteria]"),
        //         format!("        xatol         : {:?}", self.criteria.xatol),
        //         format!("        xrtol         : {:?}", self.criteria.xrtol),
        //         format!("        fatol         : {:?}", self.criteria.fatol),
        //         format!("        frtol         : {:?}", self.criteria.frtol),
        //         format!("        fltol         : {:?}", self.criteria.fltol),
        //         format!("        maxiter       : {:?}", self.criteria.maxiter),
        //         format!("        x_metric_type : {:?}", self.criteria.x_metric_type),
        //         format!("        f_metric_type : {:?}", self.criteria.f_metric_type),
        //         format!("    [Error]"),
        //         format!("        err_xa : {:?}", self.err_xa),
        //         format!("        err_xr : {:?}", self.err_xr),
        //         format!("        err_fa : {:?}", self.err_fa),
        //         format!("        err_fr : {:?}", self.err_fr),
        //         format!("        err_fl : {:?}", self.err_fl),
        //     ]
        //     .join("\n")
        // ))
        Ok(())
    }
}

// Evaluation Function
impl<X, F, J, H, M> OptimizeResult<X, F, J, H, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    J: IntoMetric<M>,
    H: IntoMetric<M>,
    M: Metric,
{
    /// Builder Pattern for setting the criteria
    pub fn set_criteria(mut self, criteria: OptimizeCriteria<X, F, M>) -> Self {
        self.criteria = criteria;
        self
    }

    /// Builder Pattern for setting the `target_f`
    pub fn set_target_f(mut self, target_f: Option<F>) -> Self {
        self.target_f = target_f;
        self
    }

    /// Increment nfev
    pub fn fev(&mut self) {
        self.nfev += 1;
    }
    /// Increment njev
    pub fn jev(&mut self) {
        self.njev += 1;
    }
    /// Increment nhev
    pub fn hev(&mut self) {
        self.nhev += 1;
    }

    /// Update new solution from solver, increment iter, and calculate error
    pub fn update(&mut self, new_solution: (X, F, Option<J>, Option<H>)) {
        (self.old_x, self.old_f) = (self.sol_x.clone(), self.sol_f.clone());
        self.sol_x = Some(new_solution.0);
        self.sol_f = Some(new_solution.1);
        self.sol_j = new_solution.2;
        self.sol_h = new_solution.3;
        self.iter += 1;
        self.err_xa = self.error_xa();
        self.err_xr = self.error_xr();
        self.err_fa = self.error_fa();
        self.err_fr = self.error_fr();
        self.err_fl = self.error_fl();
    }

    /// calculate and return `err_xa`, only calculate, without update
    fn error_xa(&self) -> Option<M> {
        Some((self.sol_x.clone()? - self.old_x.clone()?).eval(&self.criteria.x_metric_type))
    }
    /// calculate and return `err_xr`, only calculate, without update
    fn error_xr(&self) -> Option<M> {
        Some(self.error_xa()? / self.sol_x.clone()?.eval(&self.criteria.x_metric_type))
    }
    /// calculate and return `err_fa`, only calculate, without update
    fn error_fa(&self) -> Option<M> {
        Some((self.sol_f.clone()? - self.old_f.clone()?).eval(&self.criteria.f_metric_type))
    }
    /// calculate and return `err_fr`, only calculate, without update
    fn error_fr(&self) -> Option<M> {
        Some(self.error_fa()? / self.sol_f.clone()?.eval(&self.criteria.f_metric_type))
    }
    /// calculate and return `err_fl`, only calculate, without update
    fn error_fl(&self) -> Option<M> {
        Some((self.sol_f.clone()? - self.target_f.clone()?).eval(&self.criteria.f_metric_type))
    }

    /// return `true` if `xatol` in criteria is set, and satisfy
    pub fn satisfy_xatol(&self) -> bool {
        match (self.criteria.xatol, self.err_xa) {
            (Some(tol), Some(error)) => error < tol,
            _ => false,
        }
    }
    /// return `true` if `xrtol` in criteria is set, and satisfy
    pub fn satisfy_xrtol(&self) -> bool {
        match (self.criteria.xrtol, self.err_xr) {
            (Some(tol), Some(error)) => error < tol,
            _ => false,
        }
    }
    /// return `true` if `fatol` in criteria is set, and satisfy
    pub fn satisfy_fatol(&self) -> bool {
        match (self.criteria.fatol, self.err_fa) {
            (Some(tol), Some(error)) => error < tol,
            _ => false,
        }
    }
    /// return `true` if `frtol` in criteria is set, and satisfy
    pub fn satisfy_frtol(&self) -> bool {
        match (self.criteria.frtol, self.err_fr) {
            (Some(tol), Some(error)) => error < tol,
            _ => false,
        }
    }
    /// return `true` if `fltol` in criteria is set, and satisfy
    pub fn satisfy_fltol(&self) -> bool {
        match (self.criteria.fltol, self.err_fl) {
            (Some(tol), Some(error)) => error < tol,
            _ => false,
        }
    }

    /// either of all tolerance is satisfied
    pub fn satisfy_either(&self) -> bool {
        self.satisfy_xatol()
            || self.satisfy_xrtol()
            || self.satisfy_fatol()
            || self.satisfy_frtol()
            || self.satisfy_fltol()
    }

    /// return a `Vec<String>` of all the torlerance that passed
    pub fn satisfied(&self) -> Vec<String> {
        vec![
            (self.satisfy_xatol(), "xatol"),
            (self.satisfy_xrtol(), "xrtol"),
            (self.satisfy_fatol(), "fatol"),
            (self.satisfy_frtol(), "frtol"),
            (self.satisfy_fltol(), "fltol"),
        ]
        .iter()
        .filter_map(|(sat, tol)| if !sat { None } else { Some(tol.to_string()) })
        .collect::<Vec<String>>()
    }

    /// return `true` if `iter >= maxiter`, given maxiter is set in criteria
    pub fn overran(&self) -> bool {
        if let Some(maxiter) = self.criteria.maxiter {
            self.iter >= maxiter
        } else {
            false
        }
    }
    /// set flag to successful with a message
    pub fn set_success(&mut self, msg: String) {
        self.flag = OptimizeResultFlag::Success(msg);
    }
    /// set flag to failure with a message
    pub fn set_failure(&mut self, msg: String) {
        self.flag = OptimizeResultFlag::Failure(msg);
    }
}
