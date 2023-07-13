use ndarray::Array1;
use num::{complex::Complex64, Complex};
mod ba;
mod sos;
mod zpk;

#[derive(Debug, Clone, Copy)]
pub enum DesiredFilterOutput {
    Zpk,
    Ba,
    Sos,
}
/// Enum containing the filter output
#[derive(Debug, Clone)]
pub enum FilterOutput {
    /// See [Zpk]
    Zpk(Zpk),
    /// See [Ba]
    Ba(Ba),
    /// See [Sos]
    Sos(Sos),
}

/// # Zeros and poles representation
///
/// The zpk format is a 3-tuple (z, p, k), where z is an M-length array of the complex zeros of the transfer
/// function <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>z</mi><mo>=</mo><mo stretchy="false">[</mo><msub><mi>z</mi><mn>0</mn></msub><mo>,</mo><msub><mi>z</mi><mn>1</mn></msub><mo>,</mo><mo>.</mo><mo>.</mo><mo>.</mo><mo>,</mo><msub><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mi>M</mi><mo>&#x2212;</mo><mn>1</mn></mrow></msub><mo stretchy="false">]</mo></math>, p is an N-length array of the complex poles of the transfer function
/// <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>p</mi><mo>=</mo><mo stretchy="false">[</mo><msub><mi>p</mi><mn>0</mn></msub><mo>,</mo><msub><mi>p</mi><mn>1</mn></msub><mo>,</mo><mo>.</mo><mo>.</mo><mo>.</mo><mo>,</mo><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mi>N</mi><mo>&#x2212;</mo><mn>1</mn></mrow></msub><mo stretchy="false">]</mo></math>,
/// and k is a scalar gain. These represent the digital transfer function:
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>H</mi><mo stretchy="false">(</mo><mi>z</mi><mo stretchy="false">)</mo><mo>=</mo><mi>k</mi><mo>&#x22C5;</mo><mfrac><mrow><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>z</mi><mn>0</mn></msub><mo stretchy="false">)</mo><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>z</mi><mn>1</mn></msub><mo stretchy="false">)</mo><mo>&#x22EF;</mo><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>M</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msub><mo stretchy="false">)</mo></mrow><mrow><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>p</mi><mn>0</mn></msub><mo stretchy="false">)</mo><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>p</mi><mn>1</mn></msub><mo stretchy="false">)</mo><mo>&#x22EF;</mo><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>N</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msub><mo stretchy="false">)</mo></mrow></mfrac><mo>=</mo><mi>k</mi><mfrac><mrow><munderover><mo>&#x220F;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mi>M</mi><mo>&#x2212;</mo><mn>1</mn></mrow></munderover><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>z</mi><mi>i</mi></msub><mo stretchy="false">)</mo></mrow><mrow><munderover><mo>&#x220F;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mi>N</mi><mo>&#x2212;</mo><mn>1</mn></mrow></munderover><mo stretchy="false">(</mo><mi>z</mi><mo>&#x2212;</mo><msub><mi>p</mi><mi>i</mi></msub><mo stretchy="false">)</mo></mrow></mfrac></math>
///
/// or the analog transfer function:
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>H</mi><mo stretchy="false">(</mo><mi>s</mi><mo stretchy="false">)</mo><mo>=</mo><mi>k</mi><mo>&#x22C5;</mo><mfrac><mrow><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>z</mi><mn>0</mn></msub><mo stretchy="false">)</mo><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>z</mi><mn>1</mn></msub><mo stretchy="false">)</mo><mo>&#x22EF;</mo><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>M</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msub><mo stretchy="false">)</mo></mrow><mrow><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>p</mi><mn>0</mn></msub><mo stretchy="false">)</mo><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>p</mi><mn>1</mn></msub><mo stretchy="false">)</mo><mo>&#x22EF;</mo><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>p</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>N</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msub><mo stretchy="false">)</mo></mrow></mfrac><mo>=</mo><mi>k</mi><mfrac><mrow><munderover><mo>&#x220F;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mi>M</mi><mo>&#x2212;</mo><mn>1</mn></mrow></munderover><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>z</mi><mi>i</mi></msub><mo stretchy="false">)</mo></mrow><mrow><munderover><mo>&#x220F;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mrow class="MJX-TeXAtom-ORD"><mi>N</mi><mo>&#x2212;</mo><mn>1</mn></mrow></munderover><mo stretchy="false">(</mo><mi>s</mi><mo>&#x2212;</mo><msub><mi>p</mi><mi>i</mi></msub><mo stretchy="false">)</mo></mrow></mfrac><mo>.</mo></math>
///
///
/// Although the sets of roots are stored as vecs, their ordering does not matter: ([-1, -2], [-3, -4], 1) is the same filter as ([-2, -1], [-4, -3], 1).
#[derive(Debug, Clone)]
pub struct Zpk {
    pub z: Array1<Complex64>,
    pub p: Array1<Complex64>,
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

/// # Transfer function representation
///
/// The ba or tf format is a 2-tuple (b, a) representing a transfer function, where b is a length M+1 array of
/// coefficients of the M-order numerator polynomial, and a is a length N+1 array of coefficients of the N-order
/// denominator, as positive, descending powers of the transfer function variable. So the tuple of
/// <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>b</mi><mo>=</mo><mo stretchy="false">[</mo><msub><mi>b</mi><mn>0</mn></msub><mo>,</mo><msub><mi>b</mi><mn>1</mn></msub><mo>,</mo><mo>.</mo><mo>.</mo><mo>.</mo><mo>,</mo><msub><mi>b</mi><mi>M</mi></msub><mo stretchy="false">]</mo></math> and <math xmlns="http://www.w3.org/1998/Math/MathML"><mi>a</mi><mo>=</mo><mo stretchy="false">[</mo><msub><mi>a</mi><mn>0</mn></msub><mo>,</mo><msub><mi>a</mi><mn>1</mn></msub><mo>,</mo><mo>.</mo><mo>.</mo><mo>.</mo><mo>,</mo><msub><mi>a</mi><mi>N</mi></msub><mo stretchy="false">]</mo></math> can represent an analog filter of the form:
///
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>H</mi><mo stretchy="false">(</mo><mi>s</mi><mo stretchy="false">)</mo><mo>=</mo><mfrac><mrow><msub><mi>b</mi><mn>0</mn></msub><msup><mi>s</mi><mi>M</mi></msup><mo>+</mo><msub><mi>b</mi><mn>1</mn></msub><msup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>M</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msup><mo>+</mo><mo>&#x22EF;</mo><mo>+</mo><msub><mi>b</mi><mi>M</mi></msub></mrow><mrow><msub><mi>a</mi><mn>0</mn></msub><msup><mi>s</mi><mi>N</mi></msup><mo>+</mo><msub><mi>a</mi><mn>1</mn></msub><msup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>N</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msup><mo>+</mo><mo>&#x22EF;</mo><mo>+</mo><msub><mi>a</mi><mi>N</mi></msub></mrow></mfrac><mo>=</mo><mfrac><mrow><munderover><mo>&#x2211;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mi>M</mi></munderover><msub><mi>b</mi><mi>i</mi></msub><msup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>M</mi><mo>&#x2212;</mo><mi>i</mi><mo stretchy="false">)</mo></mrow></msup></mrow><mrow><munderover><mo>&#x2211;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mi>N</mi></munderover><msub><mi>a</mi><mi>i</mi></msub><msup><mi>s</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>N</mi><mo>&#x2212;</mo><mi>i</mi><mo stretchy="false">)</mo></mrow></msup></mrow></mfrac></math>
///
/// or a discrete-time filter of the form:
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>H</mi><mo stretchy="false">(</mo><mi>z</mi><mo stretchy="false">)</mo><mo>=</mo><mfrac><mrow><msub><mi>b</mi><mn>0</mn></msub><msup><mi>z</mi><mi>M</mi></msup><mo>+</mo><msub><mi>b</mi><mn>1</mn></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>M</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msup><mo>+</mo><mo>&#x22EF;</mo><mo>+</mo><msub><mi>b</mi><mi>M</mi></msub></mrow><mrow><msub><mi>a</mi><mn>0</mn></msub><msup><mi>z</mi><mi>N</mi></msup><mo>+</mo><msub><mi>a</mi><mn>1</mn></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>N</mi><mo>&#x2212;</mo><mn>1</mn><mo stretchy="false">)</mo></mrow></msup><mo>+</mo><mo>&#x22EF;</mo><mo>+</mo><msub><mi>a</mi><mi>N</mi></msub></mrow></mfrac><mo>=</mo><mfrac><mrow><munderover><mo>&#x2211;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mi>M</mi></munderover><msub><mi>b</mi><mi>i</mi></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>M</mi><mo>&#x2212;</mo><mi>i</mi><mo stretchy="false">)</mo></mrow></msup></mrow><mrow><munderover><mo>&#x2211;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mi>N</mi></munderover><msub><mi>a</mi><mi>i</mi></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo stretchy="false">(</mo><mi>N</mi><mo>&#x2212;</mo><mi>i</mi><mo stretchy="false">)</mo></mrow></msup></mrow></mfrac><mo>.</mo></math>
///
/// This “positive powers” form is found more commonly in controls engineering. If M and N are equal
/// (which is true for all filters generated by the bilinear transform), then this happens to be equivalent
/// to the “negative powers” discrete-time form preferred in DSP:
///
/// <math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mi>H</mi><mo stretchy="false">(</mo><mi>z</mi><mo stretchy="false">)</mo><mo>=</mo><mfrac><mrow><msub><mi>b</mi><mn>0</mn></msub><mo>+</mo><msub><mi>b</mi><mn>1</mn></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo>&#x2212;</mo><mn>1</mn></mrow></msup><mo>+</mo><mo>&#x22EF;</mo><mo>+</mo><msub><mi>b</mi><mi>M</mi></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo>&#x2212;</mo><mi>M</mi></mrow></msup></mrow><mrow><msub><mi>a</mi><mn>0</mn></msub><mo>+</mo><msub><mi>a</mi><mn>1</mn></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo>&#x2212;</mo><mn>1</mn></mrow></msup><mo>+</mo><mo>&#x22EF;</mo><mo>+</mo><msub><mi>a</mi><mi>N</mi></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo>&#x2212;</mo><mi>N</mi></mrow></msup></mrow></mfrac><mo>=</mo><mfrac><mrow><munderover><mo>&#x2211;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mi>M</mi></munderover><msub><mi>b</mi><mi>i</mi></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo>&#x2212;</mo><mi>i</mi></mrow></msup></mrow><mrow><munderover><mo>&#x2211;</mo><mrow class="MJX-TeXAtom-ORD"><mi>i</mi><mo>=</mo><mn>0</mn></mrow><mi>N</mi></munderover><msub><mi>a</mi><mi>i</mi></msub><msup><mi>z</mi><mrow class="MJX-TeXAtom-ORD"><mo>&#x2212;</mo><mi>i</mi></mrow></msup></mrow></mfrac><mo>.</mo></math>
///
/// Although this is true for common filters, remember that this is not true in the general case.
/// If M and N are not equal, the discrete-time transfer function coefficients must first be converted
/// to the “positive powers” form before finding the poles and zeros.
///
/// This representation suffers from numerical error at higher orders, so other formats are preferred when possible.
#[derive(Debug, Clone)]
pub struct Ba {
    pub a: Array1<Complex<f64>>,
    pub b: Array1<Complex<f64>>,
}

#[derive(Debug, Clone)]
pub struct Sos {}
