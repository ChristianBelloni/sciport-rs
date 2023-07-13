# sciport-rs

## Sciport-rs

Sciport is a collection of mathematical algorithms ported from the popular python package Scipy

## Api design

The main philosophy behind sciport is to change the api surface of scipy to better utilize the
rich rust typesystem, when deciding between keeping the original function signature and
rewriting it to better represent the valid input space, more often than not we'll decide to
change it.<br/>
for example this is the scipy butter filter api:

```python
scipy.signal.butter(N: int, Wn: array_like, btype: String, analog: bool, output: String, fs:
float)
```

Wn represents a single or a pair of frequencies and btype is the type of filter,
however, a single frequency makes sense only for a subset of btypes and so does a pair,
in our implementation we rewrite this function like:

```rust
fn filter<T>(order: u32, band_filter: BandFilter, analog: Analog) { .. }
```

where T represents the output representation of the filter (Zpk, Ba, Sos), band_filter
encapsulates the original Wn and btype like this:

```rust
enum BandFilter

pub enum BandFilter {
    Highpass(f64),
    Lowpass(f64),
    Bandpass { low: f64, high: f64 },
    Bandstop { low: f64, high: f64 },
}
```

and analog encapsulates analog and fs (since a sampling rate makes sense only when talking
about a digital filter) like this:

```rust
pub enum Analog {
    True,
    False {
        fs: f64
    }
}
```

