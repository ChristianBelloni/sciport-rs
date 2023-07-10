use std::marker::PhantomData;

use num::{Complex};



use super::{
    band_filter::BandFilter,
    iir_filter,
    output_type::{Ba, DesiredFilterOutput, FilterOutput, Zpk},
    Analog,
};

pub(crate) fn butter_filter(
    order: u32,
    band_filter: BandFilter,
    analog: Analog,
    desired_output: DesiredFilterOutput,
) -> FilterOutput {
    let proto = buttap(order);
    iir_filter(proto, order, band_filter, analog, desired_output)
}

pub struct ButterFilter<T>(PhantomData<T>);

impl ButterFilter<Zpk> {
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog) -> Zpk {
        let filter = butter_filter(order, band_filter, analog, DesiredFilterOutput::Zpk);
        filter.zpk()
    }
}

impl ButterFilter<Ba> {
    pub fn filter(order: u32, band_filter: BandFilter, analog: Analog) -> Ba {
        let filter = butter_filter(order, band_filter, analog, DesiredFilterOutput::Ba);
        filter.ba()
    }
}

fn buttap(order: u32) -> Zpk {
    let order = order as i32;
    let z = vec![];
    let range = ((-order + 1)..order).step_by(2);

    fn make_iteration(item: i32, order: i32) -> Complex<f64> {
        let numerator = Complex::new(0.0, 1.0) * std::f64::consts::PI * (item as f64);
        let denominator = 2.0 * order as f64;

        let temp = numerator / denominator;

        -temp.exp()
    }

    let p = range
        .map(|item| make_iteration(item, order))
        .collect::<Vec<_>>();

    let k = 1.0;

    Zpk { z, p, k }
}

#[cfg(test)]
mod tests {
    use super::buttap;
    use super::ButterFilter;
    use crate::signal::band_filter::BandFilter;
    use crate::signal::output_type::Ba;

    #[test]
    fn test_buttap() {
        let res = buttap(8);
        println!("{:#?}", res.p);
        assert_eq!(res.z, vec![]);
    }

    #[test]
    fn test_base_case() {
        let filter = ButterFilter::<Ba>::filter(
            8,
            BandFilter::Bandstop {
                low: 0.2,
                high: 0.4,
            },
            crate::signal::Analog::False { fs: 2.0 },
        );

        println!("{filter:#?}");
    }
}
