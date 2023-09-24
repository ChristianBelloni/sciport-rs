use std::fmt::{Debug, Display};

use num::{Float, complex::ComplexFloat};

pub mod bracket;
pub mod bisect;
pub mod regula_falsi;

pub use bracket::*;

/// The rtol for default solver
const DFAULT_RTOL    : f64 = 1e-100;

/// The ftol for default solver
const DFAULT_FTOL    : f64 = 1e-12;

/// The xtol for default solver
const DFAULT_XTOL    : f64 = 1e-12;

/// The maxiter for default solver
const DFAULT_MAXITER : u64 = 1000_u64;


/// # `RootScalarSolver<T>`
/// 
/// ## `where T: PartialEq + PartialOrder + Display + Debug`
/// The struct that hold all the terminate condition for root finding
/// 
/// 
/// ## field
/// ### `rtol: Option<T>`
/// the relative tolorence, if is `Some(_)`, the root finding algorithm will terminate if
/// `|x - x'| < x' * rtol`
/// 
/// ### `ftol: Option<T>`
/// the evaluated tolorence, if is `Some(_)`, the root finding algorithm will terminate if
/// `|f(x)| < ftol`
/// 
/// ### `xtol: Option<T>`
/// the absolute tolorence, if is `Some(_)`, the root finding algorithm will terminate if
/// `|x - x'| < xtol`
/// 
/// ### `iter: Option<u64>`
/// the maximum iteration allowed, if is `Some(_)`, the root finding algorithm will terminate if
/// `no_of_iter > maxiter`
/// 
/// 
/// ## Example
/// ### Creation
/// `RootScalarSolver<T>` implement `Default`,
/// ```
/// let solver = RootScalarSolver<f64>::default();
/// ```
/// It also utilize builder pattern which allow you to customize with ease
/// ```
/// let solver = RootScalarSolver<f64>::default()
///     .set_rtol(None)
///     .set_ftol(Some(1e-12))
///     .set_xtol(None)
///     .set_maxiter(1e10);
/// ```
/// ### Root Finding with Bracket
/// To find root for a function with bracket, you need to provide:
/// 1. A function `fun: Fn(T)->T` to find root
/// 2. the bracket `(a:T,b:T)` that the root might exist in, given that:
///     - `a` and `b` itself must be finite
///     - the function `fun` must evaluated to finite at `a` and `b`
///     - the function `fun` must evaluated to different sign at `a` and `b`
/// 3. the solving algorithm, which is an enum `BracketMethod` e.g.:
/// ```
///     BracketMethod::Bisect
///     BracketMethod::RegulaFalsi
///     BracketMethod::Ridder
///     ...
/// ````
/// #### example
/// ```
/// let solver = RootScalarSolver::<f64>::default()
///     .set_rtol(None)
///     .set_ftol(Some(1e-12));
/// 
/// fn fun(x:f64)->f64{ (x.powi(2) + 1.0_f64) * (x - 0.5_f64) }
/// let bracket = (-10_f64,10_f64);
/// let method = BracketMethod::Ridder;
/// 
/// let res: RootScalarResult<f64> = solver.solve_with_bracket(
///     &fun, bracket, method
/// );
/// 
/// let res = solver.bisect(&fun,bracket);
/// ```
pub struct RootScalarSolver<T>
where 
    T: PartialEq + PartialOrd + Display + Debug,
{
    rtol: Option<T>,
    ftol: Option<T>,
    xtol: Option<T>,
    maxiter: Option<u64>,
}

impl<T> RootScalarSolver<T>
where 
    T: PartialEq + PartialOrd + Display + Debug + Clone,
{
    pub fn set_rtol(mut self, rtol:Option<T>)->Self { self.rtol = rtol; self}
    pub fn set_ftol(mut self, ftol:Option<T>)->Self { self.ftol = ftol; self}
    pub fn set_xtol(mut self, xtol:Option<T>)->Self { self.xtol = xtol; self} 
    pub fn set_maxiter(mut self, maxiter:Option<u64>)->Self { self.maxiter = maxiter; self}
    pub fn get_rtol(&self)->Option<T> { self.rtol.clone() }
    pub fn get_ftol(&self)->Option<T> { self.ftol.clone() }
    pub fn get_xtol(&self)->Option<T> { self.xtol.clone() } 
    pub fn get_maxiter(&self)->Option<u64> { self.maxiter }
}

impl<T> Default for RootScalarSolver<T>
where 
    T: PartialEq + PartialOrd + Display + Debug + Float,
{
    fn default() -> Self{
        Self { 
            rtol: Some(T::from(DFAULT_RTOL).unwrap()), 
            ftol: Some(T::from(DFAULT_FTOL).unwrap()), 
            xtol: Some(T::from(DFAULT_XTOL).unwrap()), 
            maxiter: Some(DFAULT_MAXITER) 
        }
    }
}

impl<T> RootScalarSolver<T>
where 
    T: PartialEq + PartialOrd + Display + Debug + Float + 'static,
{

    fn validate_bracket(fun:&dyn Fn(T)->T, bracket:(T,T), res: &mut RootScalarResult<T>){        
        let (a,b) = bracket;
        if a.is_infinite() || b.is_infinite(){
            res.flag = RootScalarResultFlag::Fail("Bracket must be finite.".to_string());
            return;
        }
        let (u,v)=(fun(a),fun(b));
        if u.is_infinite() || v.is_infinite(){
            res.flag = RootScalarResultFlag::Fail("Bracket must evaluated to be finite.".to_string());
            return;
        }
        if u * v > T::zero(){
            res.flag = RootScalarResultFlag::Fail("Bracket must evaluate to different sign.".to_string());
        }
    }
    /// ### Root Finding with Bracket
    /// To find root for a function with bracket, you need to provide:
    /// 1. A function `fun: Fn(T)->T` to find root
    /// 2. the bracket `(a:T,b:T)` that the root might exist in, given that:
    ///     - `a` and `b` itself must be finite
    ///     - the function `fun` must evaluated to finite at `a` and `b`
    ///     - the function `fun` must evaluated to different sign at `a` and `b`
    ///        
    /// if any of above condition is not met, it will return result with a fail flag
    /// 
    /// 3. the solving algorithm, which is an enum `BracketMethod` e.g.:
    /// ```
    ///     BracketMethod::Bisect
    ///     BracketMethod::RegulaFalsi
    ///     BracketMethod::Ridder
    ///     ...
    /// ````
    pub fn solve_with_bracket(&self, fun:&dyn Fn(T)->T, bracket:(T,T), method: BracketMethod)->RootScalarResult<T> {
        let mut res = RootScalarResult::<T>::new();

        Self::validate_bracket(&fun, bracket.clone(), &mut res);

        let mut fun = |x|{
            res.function_called += 1;
            fun(x)
        };
        
        if res.flag.fail(){
            return res;
        }

        let mut solver = method.get_solver(&mut fun,bracket);
        let mut last_guess = None;
        
        loop{
            let (new_guess,new_y) = solver.new_guess(&mut fun);
            
            res.root = new_guess;
            res.error = new_y;

            if new_y.is_infinite(){
                res.flag = RootScalarResultFlag::Fail("Function evaluated to be infinite".to_string());
                break
            }

            if let Some(ftol) = self.ftol{
                if new_y.abs() < ftol {
                    res.flag.ftol_reach();
                }
            }
            if let Some(last_guess) = last_guess{
                if let Some(rtol) = self.rtol{
                    if (new_guess-last_guess).abs() < new_guess * rtol{
                        res.flag.rtol_reach();
                    }
                }
    
                if let Some(xtol) = self.xtol{
                    if (last_guess-new_guess).abs() < xtol{
                        res.flag.xtol_reach();
                    }
                }
            }
            
            if res.flag.sucess(){break;}

            if let Some(maxiter) = self.maxiter{
                if res.iterations >= maxiter{
                    res.flag = RootScalarResultFlag::Fail("Fail, max iterations reached".to_string());
                    break;
                }
            }
            res.iterations += 1;
            last_guess = Some(new_guess)
        }
        res
    }

    pub fn bisect(&self,fun:&dyn Fn(T)->T,bracket:(T,T))->RootScalarResult<T>{
        self.solve_with_bracket(fun, bracket, BracketMethod::Bisect)
    }
    pub fn regula_falsi(&self,fun:&dyn Fn(T)->T,bracket:(T,T))->RootScalarResult<T>{
        self.solve_with_bracket(fun, bracket, BracketMethod::RegulaFalsi)
    }
    pub fn ridder(&self,fun:&dyn Fn(T)->T,bracket:(T,T))->RootScalarResult<T>{
        self.solve_with_bracket(fun, bracket, BracketMethod::Ridder)
    }

}


impl<T> RootScalarSolver<T>
where 
    T: PartialEq + PartialOrd + Display + Debug + ComplexFloat,
{}


/// ## `RootScalarResult` flag
/// ```
/// Running
/// ```
/// running, should not be seen, just for internal logic
/// 
/// ```
/// Fail(String)
/// ```
/// root finding fail, and reason is stated in the `String`
/// 
/// ```
/// Success{rot:bool,ftol:bool,xtol:bool}
/// ```
/// root finding sucess, struct body specify which condition is reached
#[derive(Debug,Clone)]
pub enum RootScalarResultFlag{
    Running,
    Success{
        rtol: bool,
        ftol: bool,
        xtol: bool,
    },
    Fail(String),
}

impl RootScalarResultFlag{
    /// return `true` if the flag is `Sucess{_}`
    pub fn sucess(&self)->bool{
        match self{
            Self::Success{rtol:_,ftol:_,xtol:_}=>true,
            _=>false
        }
    }
    /// return `true` if the flag is `Fail(String)`
    pub fn fail(&self)->bool{
        match self{
            Self::Fail(_)=>true,
            _=>false
        }
    }
    /// Set `rtol` in `Success` flag to be true,
    /// if the flag is no `Success` it change it to `Success`
    pub fn rtol_reach(&mut self){
        match self {
            Self::Success { rtol:_, ftol, xtol }=>{
                *self=Self::Success { rtol: true, ftol:*ftol, xtol:*xtol};
            }
            _=>{
                *self=Self::Success { rtol: true, ftol: false, xtol: false };
            }
        }
    }
    /// Set `ftol` in `Success` flag to be true,
    /// if the flag is no `Success` it change it to `Success`
    pub fn ftol_reach(&mut self){
        match self {
            Self::Success { rtol, ftol:_, xtol }=>{
                *self=Self::Success { rtol: *rtol, ftol:true, xtol:*xtol};
            }
            _=>{
                *self=Self::Success { rtol: false, ftol: true, xtol: false };
            }
        }
    }
    /// Set `xtol` in `Success` flag to be true,
    /// if the flag is no `Success` it change it to `Success`
    pub fn xtol_reach(&mut self){
        match self {
            Self::Success { rtol, ftol, xtol:_ }=>{
                *self=Self::Success { rtol: *rtol, ftol:*ftol, xtol:true};
            }
            _=>{
                *self=Self::Success { rtol: false, ftol: false, xtol: true };
            }
        }
    }
}

/// # `RootScalarResult<T>`
/// ## where T: PartialEq + PartialOrd + Display + Debug
/// Store the result of root finding
/// ## field
/// ### `root: T`
/// the result, the last guess of the algorithm, 
/// check the flag to see if the algorithm success
/// 
/// ### `error: T`
/// the evaluated error of the guess.
/// 
/// ### `iterations: u64`
/// number of iterations executed.
/// 
/// ### `function_called: u64`
/// number of function called.
/// 
/// ### `flag: RootScalarResultFlag`
/// the flag for of the result, 
/// indicate the success and fail condition
pub struct RootScalarResult<T>
where
    T: PartialEq + PartialOrd + Display + Debug,
{
    pub root: T,
    pub error: T,
    pub iterations: u64,
    pub function_called: u64,
    pub flag: RootScalarResultFlag,
}

impl<T>  RootScalarResult<T>
where
    T: PartialEq + PartialOrd + Display + Debug + Float,
{
    fn new() -> Self {
        Self { 
            root: T::zero(), 
            error: T::zero(),
            iterations: 0, 
            function_called: 0, 
            flag: RootScalarResultFlag::Running 
        }
    }
}

impl<T> ToString for RootScalarResult<T>
where
    T: PartialEq + PartialOrd + Display + Debug,
{
    fn to_string(&self) -> String {
        vec![
            format!("root : {:>+10.4}", self.root),
            format!("error : {:>+10.4}", self.error),
            format!("iteration : {:>6}", self.iterations),
            format!("no. of function called : {:>6}", self.function_called),
            format!("flag : {:?}", self.flag),
        ].join(", ")
    }
}

impl<T> Debug for RootScalarResult<T> 
where
    T: PartialEq + PartialOrd + Display + Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"{}", self.to_string())
    }    
}