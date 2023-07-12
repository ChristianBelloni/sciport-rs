use num::Complex;
mod ba;
mod sos;
mod zpk;

#[derive(Debug, Clone, Copy)]
pub enum DesiredFilterOutput {
    Zpk,
    Ba,
    Sos,
}

#[derive(Debug, Clone)]
pub enum FilterOutput {
    Zpk(Zpk),
    Ba(Ba),
    Sos(Sos),
}

#[derive(Debug, Clone)]
pub struct Zpk {
    pub z: Vec<Complex<f64>>,
    pub p: Vec<Complex<f64>>,
    pub k: f64,
}

impl FilterOutput {
    pub fn new(data: Zpk) -> Self {
        Self::Zpk(data)
    }

    pub fn zpk(self) -> Zpk {
        match self {
            Self::Zpk(data) => data,
            _ => unreachable!(),
        }
    }

    pub fn ba(self) -> Ba {
        match self {
            Self::Ba(data) => data,
            _ => unreachable!(),
        }
    }

    pub fn sos(self) -> Sos {
        match self {
            Self::Sos(data) => data,
            _ => unreachable!(),
        }
    }

    pub fn get_output(input: Zpk, desired: DesiredFilterOutput) -> Self {
        match desired {
            DesiredFilterOutput::Zpk => Self::Zpk(input),
            DesiredFilterOutput::Ba => Self::Ba(input.into()),
            _ => todo!(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Ba {
    pub a: Vec<Complex<f64>>,
    pub b: Vec<Complex<f64>>,
}

#[derive(Debug, Clone)]
pub struct Sos {}
