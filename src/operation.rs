//! Per-bin combination: `Operation` enum and `getRatio` kernel.

use crate::CompareOpts;

/// How two scaled per-bin coverage values are combined.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    Log2,
    Ratio,
    ReciprocalRatio,
    Subtract,
    Add,
    Mean,
    First,
    Second,
}

impl std::str::FromStr for Operation {
    type Err = String;
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "log2" => Ok(Self::Log2),
            "ratio" => Ok(Self::Ratio),
            "reciprocal_ratio" => Ok(Self::ReciprocalRatio),
            "subtract" => Ok(Self::Subtract),
            "add" => Ok(Self::Add),
            "mean" => Ok(Self::Mean),
            "first" => Ok(Self::First),
            "second" => Ok(Self::Second),
            _ => Err(format!(
                "unknown operation '{s}'; choose log2 ratio reciprocal_ratio \
                 subtract add mean first second"
            )),
        }
    }
}

impl Operation {
    /// Whether this op is ratio-family (pseudocount applies only to these).
    pub(crate) fn is_ratio(self) -> bool {
        matches!(self, Self::Log2 | Self::Ratio | Self::ReciprocalRatio)
    }
}

/// deeptools `getRatio`: scale, NaN-propagate, then combine.
pub(crate) fn get_ratio(cov1: f64, cov2: f64, opts: &CompareOpts) -> f64 {
    let v1 = opts.scale_factors[0] * cov1;
    let v2 = opts.scale_factors[1] * cov2;
    if v1.is_nan() || v2.is_nan() {
        return f64::NAN;
    }
    if opts.operation.is_ratio() {
        let num = v1 + opts.pseudocount[0];
        let den = v2 + opts.pseudocount[1];
        let ratio = num / den;
        match opts.operation {
            Operation::Log2 => ratio.log2(),
            Operation::Ratio => ratio,
            Operation::ReciprocalRatio => {
                if ratio >= 1.0 {
                    ratio
                } else {
                    -1.0 / ratio
                }
            }
            _ => unreachable!(),
        }
    } else {
        match opts.operation {
            Operation::Subtract => v1 - v2,
            Operation::Add => v1 + v2,
            // deeptools/numpy compute `(v1 + v2) / 2.0` (not f64::midpoint);
            // the bedGraph must match it bit-for-bit.
            #[allow(clippy::manual_midpoint)]
            Operation::Mean => (v1 + v2) / 2.0,
            Operation::First => v1,
            Operation::Second => v2,
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(op: Operation) -> CompareOpts {
        CompareOpts {
            operation: op,
            ..CompareOpts::default()
        }
    }

    #[test]
    fn ratio_doctest() {
        // deeptools getRatio doctest: ratio [9,19] pc [1,1] = 10/20 = 0.5
        let mut o = opts(Operation::Ratio);
        o.pseudocount = [1.0, 1.0];
        assert!((get_ratio(9.0, 19.0, &o) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn ratio_zero_over_zero() {
        // [0,0] pc [1,1] = 1/1 = 1.0
        let o = opts(Operation::Ratio);
        assert!((get_ratio(0.0, 0.0, &o) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nan_propagates() {
        let o = opts(Operation::Ratio);
        assert!(get_ratio(f64::NAN, 1.0, &o).is_nan());
        assert!(get_ratio(1.0, f64::NAN, &o).is_nan());
    }

    #[test]
    fn subtract_no_pseudocount() {
        let o = opts(Operation::Subtract);
        assert!((get_ratio(20.0, 10.0, &o) - 10.0).abs() < 1e-12);
    }

    #[test]
    fn reciprocal_ratio_no_pc() {
        let mut o = opts(Operation::ReciprocalRatio);
        o.pseudocount = [0.0, 0.0];
        assert!((get_ratio(2.0, 1.0, &o) - 2.0).abs() < 1e-12);
        assert!((get_ratio(1.0, 2.0, &o) - (-2.0)).abs() < 1e-12);
        assert!((get_ratio(1.0, 1.0, &o) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn scale_factor_applied() {
        let mut o = opts(Operation::First);
        o.scale_factors = [0.5, 1.0];
        assert!((get_ratio(10.0, 20.0, &o) - 5.0).abs() < 1e-12);
    }
}
