use std::fmt::{Display, Debug};
use num::Float;

use crate::optimize::root_scalar::BracketSolver;

/// ## Bisect Algorithm
pub struct BisectSolver<T>{
    pos:T,neg:T
}

impl<T> BracketSolver<T> for BisectSolver<T>
where
    T: PartialEq + PartialOrd + Display + Debug + Float,
{
    fn init(fun:&mut dyn FnMut(T)->T, bracket:(T,T))->Box<Self> {
        let (pos, neg) = {
            if fun(bracket.0).is_sign_positive() {bracket}
            else{ (bracket.1,bracket.0) }
        };
        Box::new(Self { pos, neg })
    }
    fn new_guess(&mut self,fun:&mut dyn FnMut(T)->T)->(T,T){
        let mid = (self.pos+self.neg)/T::from(2.0).unwrap();
        let y = fun(mid);
        if y.is_sign_positive(){ self.pos = mid; }
        else{ self.neg = mid; }
        (mid,y)
    }
}
