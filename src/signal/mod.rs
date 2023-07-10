use self::{
    band_filter::{lp2bf_zpk, BandFilter},
    output_type::{DesiredFilterOutput, FilterOutput, Zpk},
    tools::bilinear_zpk,
};

pub mod band_filter;
pub mod bessel;
pub mod butter;
pub mod output_type;
pub mod tools;

#[derive(Debug, Clone, Copy)]
pub enum Analog {
    True,
    False { fs: f64 },
}

impl Analog {
    pub fn is_analog(&self) -> bool {
        match self {
            Self::True => true,
            Self::False { .. } => false,
        }
    }
}

pub(crate) fn iir_filter(
    proto: Zpk,
    _order: u32,
    mut band_filter: BandFilter,
    mut analog: Analog,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let mut warped = band_filter;
    match &mut analog {
        Analog::True => {}
        Analog::False { fs } => {
            band_filter = (band_filter * 2.0) / *fs;
            *fs = 2.0;
            warped = ((band_filter * std::f64::consts::PI) / *fs).tan() * 2.0 * *fs;
        }
    }

    let mut result = lp2bf_zpk(proto, warped);
    if let Analog::False { fs } = &analog {
        result = bilinear_zpk(result, *fs);
    }

    FilterOutput::get_output(result, desired_output)
}
