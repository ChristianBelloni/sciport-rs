use crate::optimize::*;

/// Criteria
#[derive(Debug, Clone)]
pub struct OptimizeCriteria<X, F, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    M: Metric,
{
    /// Satisfies xatol if `|x-x'| < xatol`
    pub xatol: Option<M>,
    /// Satisfies xrtol if `|x-x'| < xrtol * x'`
    pub xrtol: Option<M>,
    /// Satisfies fatol if `|f-f'| < fatol`
    pub fatol: Option<M>,
    /// Satisfies frtol if `|f-f'| < frtol * f'`
    pub frtol: Option<M>,
    /// Satisfies fltol if `|f-target_f| < fatol`
    pub fltol: Option<M>,
    /// Fail if `iter > maxiter`
    pub maxiter: Option<u64>,
    /// specify the metric evaluation type for x
    pub x_metric_type: MetricType<X, M>,
    /// specify the metric evaluation type for f
    pub f_metric_type: MetricType<F, M>,
}

/// Default `xatol`
const DEFAULT_XATOL: f64 = 1e-9;
/// Default `xrtol`
const DEFAULT_XRTOL: f64 = 1e-100;
/// Default `fatol`
const DEFAULT_FATOL: f64 = 1e-9;
/// Default `frtol`
const DEFAULT_FRTOL: f64 = 1e-100;
/// Default `fltol`
const DEFAULT_FLTOL: f64 = 1e-9;
/// Default `maxiter`
const DEFAULT_MAXITER: u64 = 1000;

impl<X, F, M> OptimizeCriteria<X, F, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    M: Metric,
{
    ///  Builder Pattern for setting `xatol`
    pub fn set_xatol(mut self, tol: Option<M>) -> Self {
        self.xatol = tol;
        self
    }
    ///  Builder Pattern for setting `xrtol`
    pub fn set_xrtol(mut self, tol: Option<M>) -> Self {
        self.xrtol = tol;
        self
    }
    ///  Builder Pattern for setting `fatol`
    pub fn set_fatol(mut self, tol: Option<M>) -> Self {
        self.fatol = tol;
        self
    }
    ///  Builder Pattern for setting `frtol`
    pub fn set_frtol(mut self, tol: Option<M>) -> Self {
        self.frtol = tol;
        self
    }
    ///  Builder Pattern for setting `fltol`
    pub fn set_fltol(mut self, tol: Option<M>) -> Self {
        self.fltol = tol;
        self
    }
    ///  Builder Pattern for setting `maxiter`
    pub fn set_maxiter(mut self, max: Option<u64>) -> Self {
        self.maxiter = max;
        self
    }
    ///  Builder Pattern for setting `x_metric_type`
    pub fn set_x_metric_type(mut self, metric_type: MetricType<X, M>) -> Self {
        self.x_metric_type = metric_type;
        self
    }
    ///  Builder Pattern for setting `f_metric_type`
    pub fn set_f_metric_type(mut self, metric_type: MetricType<F, M>) -> Self {
        self.f_metric_type = metric_type;
        self
    }

    /// Create a new criteria with no parameter, and default `Metric::L2Norm` for both x and f
    pub fn empty() -> Self {
        Self {
            xatol: None,
            xrtol: None,
            fatol: None,
            frtol: None,
            fltol: None,
            maxiter: None,
            x_metric_type: MetricType::L2Norm,
            f_metric_type: MetricType::L2Norm,
        }
    }
}

impl<X, F, M> Default for OptimizeCriteria<X, F, M>
where
    X: IntoMetric<M>,
    F: IntoMetric<M>,
    M: Metric,
{
    fn default() -> Self {
        OptimizeCriteria {
            xatol: M::from(DEFAULT_XATOL),
            xrtol: M::from(DEFAULT_XRTOL),
            fatol: M::from(DEFAULT_FATOL),
            frtol: M::from(DEFAULT_FRTOL),
            fltol: M::from(DEFAULT_FLTOL),
            maxiter: Some(DEFAULT_MAXITER),
            x_metric_type: MetricType::L2Norm,
            f_metric_type: MetricType::L2Norm,
        }
    }
}
