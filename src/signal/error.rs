use std::fmt::Debug;
use thiserror::Error;

use super::filter_design;

#[derive(Debug, Error)]
pub enum Error {
    #[error(transparent)]
    IIRFilter(#[from] filter_design::error::Error),
}
