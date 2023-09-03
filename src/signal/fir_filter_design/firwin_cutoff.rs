use ndarray::Array1;

#[derive(Debug, Clone)]
pub enum FirwinCutoff {
    Freq(f64),
    BandEdges(Array1<f64>),
}

impl std::ops::Div<f64> for FirwinCutoff {
    type Output = Self;
    fn div(self, rhs: f64) -> Self::Output {
        match self {
            FirwinCutoff::Freq(a) => Self::Freq(a / rhs),
            FirwinCutoff::BandEdges(a) => Self::BandEdges(a / rhs),
        }
    }
}
