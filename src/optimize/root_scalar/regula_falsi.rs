use std::fmt::{Display, Debug};
use num::Float;

use crate::optimize::root_scalar::BracketSolver;

pub struct RegulaFalsiSolver<T>{
    pos:T,
    neg:T,
    pos_y:T,
    neg_y:T,
}

impl<T> BracketSolver<T> for RegulaFalsiSolver<T>
where
    T: PartialEq + PartialOrd + Display + Debug + Float,
{
    fn init(fun:&mut dyn FnMut(T)->T, bracket:(T,T))->Box<Self> {
        let (pos, neg,pos_y,neg_y) = {
            let (x1,x2) = bracket;
            let (y1,y2) = (fun(x1),fun(x2));
            if y1.is_sign_positive() {(x1,x2,y1,y2)}
            else{ (x2,x1,y2,y1) }
        };
        Box::new(Self { pos, neg, pos_y, neg_y})
    }
    fn new_guess(&mut self,fun:&mut dyn FnMut(T)->T)->(T,T){
        let mid = (self.pos_y*self.neg-self.neg_y*self.pos)/(self.pos_y-self.neg_y+T::epsilon());
        
        let y = fun(mid);
        
        if y.is_sign_positive(){ self.pos = mid;self.pos_y = y; }
        else{ self.neg = mid; self.neg_y = y; }
        
        (mid,y)
    }
}

pub struct RidderSolver<T>{
    x0:T,
    x2:T,
    y0:T,
    y2:T,
}

impl<T> BracketSolver<T> for RidderSolver<T>
where
    T: PartialEq + PartialOrd + Display + Debug + Float,
{
    /*
    Ridders, C. (1979). 
    "A new algorithm for computing a single root of a real continuous function". 
    IEEE Transactions on Circuits and Systems. 26: 979–980. doi:10.1109/TCS.1979.1084580
    */
    fn init(fun:&mut dyn FnMut(T)->T, bracket:(T,T))->Box<Self> {
    let (x0,x2) = bracket;
        let (y0,y2) = (fun(x0),fun(x2));
        Box::new(Self { x0,x2,y0,y2})
    }

    fn new_guess(&mut self,fun:&mut dyn FnMut(T)->T)->(T,T){
        let x1 = (self.x0 + self.x2)/ T::from(2.0).unwrap();
        let d = self.x0 - x1;

        let y1 = fun(x1);

        let x3 = x1 + d * (y1/self.y2) / (((y1/self.y0).powi(2)-self.y2/self.y0).sqrt()+T::epsilon());
        let y3 = fun(x3);

        if (y1/y3).is_sign_negative(){
            self.x0 = x1;
            self.y0 = y1;
            self.x2 = x3;
            self.y2 = y3;
        }
        else if (self.y2/y3).is_sign_negative(){
            self.x0 = x3;
            self.y0 = y3;
        }
        else{
            self.x2 = x3;
            self.y2 = y3;
        }

        (x3,y3)
    }

}