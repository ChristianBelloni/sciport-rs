use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    Bessel(#[from] super::bessel::Error),
    #[error(transparent)]
    ButterOrd(#[from] super::butterord::Error),
    #[error("{0}")]
    Infallible(#[from] Infallible),
}

#[derive(Debug, Error)]
pub enum Infallible {}
