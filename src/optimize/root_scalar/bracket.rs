use num::Float;


use crate::optimize::root_scalar::*;
use crate::optimize::root_scalar::bisect::*;
use crate::optimize::root_scalar::regula_falsi::*;

/// Specify different type of root finding algorithms
#[derive(Debug,Clone, Copy)]
pub enum BracketMethod{
    /// ## Bisect Algorithm
    /// repeatly find the mid point between bracket,
    /// and progressively narrow the bracket to find root by checking the sign
    Bisect,
    /// ## Regula Falsi Algorithm/ False Position Algorithm
    /// Repeatly find the x-intercept `c` of the line that pass througth`(a,f(a))` and `(b,f(b))`
    /// and narrow the bracket by checking the sign of `f(c)`
    RegulaFalsi,
    /// ## Ridder Algorithm
    /// Base on Regula Falsi.
    /// 
    /// Ridders, C. (1979). 
    /// "A new algorithm for computing a single root of a real continuous function". 
    /// IEEE Transactions on Circuits and Systems. 26: 979–980. doi:10.1109/TCS.1979.1084580
    Ridder,
    Brent,
}

impl BracketMethod{
    pub fn get_solver<T>(&self,fun:&mut dyn FnMut(T)->T, bracket:(T,T))-> Box<dyn BracketSolver<T>>
    where
        T: PartialEq + PartialOrd + Display + Debug + Float + 'static,
    {
        match self {
            BracketMethod::Bisect=>BisectSolver::init(fun,bracket),
            BracketMethod::RegulaFalsi=>RegulaFalsiSolver::init(fun, bracket),
            BracketMethod::Ridder=>RidderSolver::init(fun, bracket),
            _=>unreachable!(),
        }
    }
}

pub trait BracketSolver<T>
where
    T: PartialEq + PartialOrd + Display + Debug + Float,
{
    fn init(fun:&mut dyn FnMut(T)->T, bracket:(T,T))->Box<Self> where Self:Sized;
    fn new_guess(&mut self,fun:&mut dyn FnMut(T)->T)->(T,T);
}

