use std::fmt::{Debug, Display};
use std::ops::{Add, Sub};
use std::rc::Rc;

use num::complex::{Complex32, Complex64, ComplexFloat};
use num::traits::FloatConst;
use num::{Complex, Float, NumCast, One, Zero};

/// # `Metric`
/// `Metric` is trait for float, which all compare the optimizing solution to the allowed tolerance
///
/// It is implemented for `f32` and `f64`
pub trait Metric: Float + Sized + Clone {}
impl Metric for f32 {}
impl Metric for f64 {}

/// # MetricType
/// Different type of method to measure the norm of a certain type
#[derive(Clone)]
pub enum MetricType<T, M>
where
    T: IntoMetric<M>,
    M: Metric,
{
    /// powered sum of all element
    PowerSum(M),
    /// L1 norm
    L1Norm,
    /// L2 norm
    L2Norm,
    /// p norm
    PNorm(M),
    /// mean square
    MS,
    /// root mean square
    RMS,
    /// custom function
    Custom(Rc<dyn Fn(&T) -> M>),
}

impl<T, M> Debug for MetricType<T, M>
where
    T: IntoMetric<M>,
    M: Metric + Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_fmt(format_args!("MetricType::{}", {
            match self {
                Self::PNorm(p) => format!("PNorm-{:?}", p),
                Self::L1Norm => format!("L1-Norm"),
                Self::L2Norm => format!("L2-Norm"),
                Self::MS => format!("MS"),
                Self::RMS => format!("RMS"),
                Self::PowerSum(p) => format!("PowerSum-{:?}", p),
                Self::Custom(_) => format!("Custom"),
            }
        }))
    }
}

/// # `IntoMetric<M>`
/// `IntoMetric<M>` is a trait for evaluating optimization solution.
///
/// Since optimization solution is not always comparable with the specified tolerance,
/// e.g. the solution with type `Array<f64>` and tolerance metric with type `f64`
///
/// this trait allow all optimizaition solution to be turn into a `Metric`. Implemented for:
/// ```ignore
/// f32, f64, Complex32, Complex64, Array1<f32>, Array1<f64>, Array1<Complex32>, Array1<Complex64>
/// ```
///
/// Type with this trait must meet bound `Sub<Output=Self>`,
/// since its nesserary to compare new and old solution in iterative optimization.
pub trait IntoMetric<M>
where
    Self: Sub<Output = Self> + Add<Output = Self> + Sized + Clone,
    M: Metric,
{
    /// return the total number of element in the type for calculating mean.
    fn n(&self) -> M;
    /// return the sum of powered `p` of element as a metric.
    fn power_sum(&self, p: M) -> M;
    /// return the L1 norm as a metric.
    fn l1_norm(&self) -> M {
        self.p_norm(M::one())
    }
    /// return the L2 norm as a metric.
    fn l2_norm(&self) -> M {
        self.p_norm(M::from(2).unwrap())
    }
    /// return the p norm as a metric.
    fn p_norm(&self, p: M) -> M {
        self.power_sum(p).powf(p.recip())
    }
    /// return the mean square as a metric.
    fn ms(&self) -> M {
        self.power_sum(M::from(2).unwrap()) / self.n()
    }
    /// return the root mean square as a metric.
    fn rms(&self) -> M {
        self.l2_norm() / self.n().powi(2)
    }
    /// return evaluate by specified metric type as a metric.
    fn eval(&self, m: &MetricType<Self, M>) -> M
    where
        Self: Sized,
    {
        match m {
            MetricType::PowerSum(p) => self.power_sum(*p),
            MetricType::L1Norm => self.l2_norm(),
            MetricType::L2Norm => self.l1_norm(),
            MetricType::PNorm(p) => self.p_norm(*p),
            MetricType::MS => self.ms(),
            MetricType::RMS => self.rms(),
            MetricType::Custom(f) => f(self),
        }
    }
}

impl<S> IntoMetric<S> for S
where
    S: Float + Sized + Clone,
    S: Metric,
{
    fn n(&self) -> S {
        S::one()
    }

    fn power_sum(&self, p: S) -> S {
        self.abs().powf(p)
    }
}

impl<S> IntoMetric<S> for Complex<S>
where
    S: Float + FloatConst + Sized + Clone,
    S: Metric,
{
    fn n(&self) -> S {
        S::one()
    }

    fn power_sum(&self, p: S) -> S {
        S::from(self.abs().powf(p)).unwrap()
    }
}

impl<S, M> IntoMetric<M> for Array1<S>
where
    S: IntoMetric<M>,
    M: Metric,
{
    fn n(&self) -> M {
        M::from(self.len()).unwrap()
    }

    fn power_sum(&self, p: M) -> M {
        self.fold(M::zero(), |prd, x| prd + x.power_sum(p))
    }
}

/// A macro for implementing `IntoMetric` for `f32`,`f64`, `Complex32` and `Complexf64`, into metric `f32` and `f64` respectively.
macro_rules! impl_metric_complexfloat {
    ($trait_name:ident,$type_name:ident, $t:ty, $m:ty) => {
        impl IntoMetric<$m> for $t {
            fn n(&self) -> $m {
                <$m>::one()
            }
            fn power_sum(&self, p: $m) -> $m {
                (self).abs().powf(p)
            }
        }
        pub type $type_name = MetricType<$t, $m>;
        pub trait $trait_name: IntoMetric<$m> {}
    };
}

// impl_metric_complexfloat!(R32IntoMetric, R32MetricType, f32, f32);
// impl_metric_complexfloat!(R64IntoMetric, R64MetricType, f64, f64);
// impl_metric_complexfloat!(Z32IntoMetric, Z32MetricType, Complex32, f32);
// impl_metric_complexfloat!(Z64IntoMetric, Z64MetricType, Complex64, f64);

use ndarray::Array1;

/// A macro for implementing `IntoMetric` for `Array1` with type `f32`,`f64`
/// , `Complex32` and `Complexf64`, into metric `f32` and `f64` respectively.
macro_rules! impl_metric_array1 {
    ($trait_name:ident,$type_name:ident, $t:ty, $m:tt) => {
        impl IntoMetric<$m> for $t {
            fn n(&self) -> $m {
                $m::from(self.len() as u8)
            }
            fn power_sum(&self, p: $m) -> $m {
                self.fold($m::zero(), |prd, x| prd + x.abs().powf(p))
            }
        }
        pub type $type_name = MetricType<$t, $m>;
        pub trait $trait_name: IntoMetric<$m> {}
    };
}

// impl_metric_array1!(ArrayR32IntoMetric, Array1R32MetricType, Array1<f32>, f32);
// impl_metric_array1!(ArrayR64IntoMetric, Array1R64MetricType, Array1<f64>, f64);
// impl_metric_array1!(
//     ArrayZ32IntoMetric,
//     Array1Z32MetricType,
//     Array1<Complex32>,
//     f32
// );
// impl_metric_array1!(
//     ArrayZ64IntoMetric,
//     Array1Z64MetricType,
//     Array1<Complex64>,
//     f64
// );
