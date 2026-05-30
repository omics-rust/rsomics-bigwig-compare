//! bedGraph emission: per-chromosome writer, run-length merge, Python `{:g}` formatter.

use std::io::Write;

use rsomics_bbi::BigWig;
use rsomics_common::{Result, RsomicsError};

use crate::CompareOpts;
use crate::bins::binned_values;
use crate::operation::get_ratio;

/// Write one chromosome's combined bins, mirroring deeptools
/// `writeBedGraph_worker`: skip-zero-over-zero, then either fixedStep (every
/// bin) or run-length merge (adjacent equal values; NaN bins never written; a
/// trailing run whose value is 0 or NaN is also not written).
pub(crate) fn write_chrom(
    out: &mut impl Write,
    a: &mut BigWig,
    b: &mut BigWig,
    chrom: &str,
    len: u32,
    opts: &CompareOpts,
) -> Result<u64> {
    let cov1 = binned_values(a, chrom, len, opts.bin_size, opts.missing_data_as_zero)?;
    let cov2 = binned_values(b, chrom, len, opts.bin_size, opts.missing_data_as_zero)?;

    let bin_size = u64::from(opts.bin_size);
    let chrom_len = u64::from(len);
    let n = cov1.len();
    let mut lines: u64 = 0;

    let mut prev: Option<f64> = None;
    let mut write_start: u64 = 0;
    let mut write_end: u64 = 0;

    for i in 0..n {
        let c1 = cov1[i];
        let c2 = cov2[i];

        // skipZeroOverZero: sum of the two coverage values == 0, before
        // pseudocount. NaN never equals 0 so a NaN bin is not skipped here.
        if opts.skip_zero_over_zero && (c1 + c2) == 0.0 {
            prev = None;
            continue;
        }

        let value = get_ratio(c1, c2, opts);

        if opts.fixed_step {
            let ws = i as u64 * bin_size;
            let we = (ws + bin_size).min(chrom_len);
            write_line(out, chrom, ws, we, value)?;
            lines += 1;
            continue;
        }

        match prev {
            None => {
                write_start = i as u64 * bin_size;
                write_end = (write_start + bin_size).min(chrom_len);
                prev = Some(value);
            }
            Some(pv) if bits_eq(pv, value) => {
                write_end = (write_end + bin_size).min(chrom_len);
            }
            Some(pv) => {
                if !pv.is_nan() {
                    write_line(out, chrom, write_start, write_end, pv)?;
                    lines += 1;
                }
                prev = Some(value);
                write_start = write_end;
                write_end = (write_start + bin_size).min(chrom_len);
            }
        }
    }

    if !opts.fixed_step
        && let Some(pv) = prev
        // deeptools: `if previousValue and writeStart != end and not isnan`.
        // `previousValue` is falsy for 0.0, so a trailing 0-run is not written.
        && pv != 0.0
        && !pv.is_nan()
        && write_start != chrom_len
    {
        write_line(out, chrom, write_start, chrom_len, pv)?;
        lines += 1;
    }

    Ok(lines)
}

fn write_line(out: &mut impl Write, chrom: &str, start: u64, end: u64, value: f64) -> Result<()> {
    let s = format_g(value);
    writeln!(out, "{chrom}\t{start}\t{end}\t{s}").map_err(RsomicsError::Io)
}

/// deeptools merges on Python `==`, which is exact float equality (NaN never
/// equal). We mirror that: bit-exact equality, NaN != NaN.
fn bits_eq(a: f64, b: f64) -> bool {
    a == b
}

/// Format a float like Python's `{:g}` (6 significant digits, trailing zeros stripped).
pub(crate) fn format_g(v: f64) -> String {
    if v == 0.0 {
        // {:g} prints negative-zero as "-0"; numpy/deeptools emit "0".
        return "0".to_owned();
    }
    if v.is_nan() {
        return "nan".to_owned();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf" } else { "-inf" }.to_owned();
    }
    python_g(v)
}

/// Python `{:g}`: 6 significant digits, switching to exponent form outside
/// `1e-4..1e16`, with trailing zeros (and a bare trailing `.`) stripped.
///
/// The decimal exponent is taken from a 6-sig-fig scientific render
/// (`{:.5e}`) rather than `log10().floor()` — that render already rounds to
/// the precision Python uses, so its exponent is the post-rounding one (a
/// value like `999999.6` rounds up to `1e6`, which `log10().floor()` would
/// mis-bucket as exponent 5).
fn python_g(v: f64) -> String {
    let sci = format!("{v:.5e}");
    let (_, exp_str) = sci.split_once('e').unwrap();
    let exp: i32 = exp_str.parse().unwrap();

    if !(-4..16).contains(&exp) {
        return normalise_exponential(&sci);
    }
    let decimals = usize::try_from((5 - exp).max(0)).unwrap();
    let s = format!("{v:.decimals$}");
    let s = if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        &s
    };
    s.to_owned()
}

/// Rust's `{:e}` gives e.g. `1.5e2` / `1e-5`; Python `{:g}` wants `1.5e+02` /
/// `1e-05` (sign always present, exponent ≥ 2 digits) with mantissa zeros
/// stripped.
fn normalise_exponential(s: &str) -> String {
    let (mantissa, exp) = s.split_once('e').unwrap();
    let mantissa = if mantissa.contains('.') {
        mantissa.trim_end_matches('0').trim_end_matches('.')
    } else {
        mantissa
    };
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(rest) => ('-', rest),
        None => ('+', exp.strip_prefix('+').unwrap_or(exp)),
    };
    format!("{mantissa}e{sign}{digits:0>2}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g_format_basics() {
        assert_eq!(format_g(0.0), "0");
        assert_eq!(format_g(-0.0), "0");
        assert_eq!(format_g(1.0), "1");
        assert_eq!(format_g(0.5), "0.5");
        assert_eq!(format_g(f64::NAN), "nan");
        assert_eq!(format_g(1.5), "1.5");
        assert_eq!(format_g(-1.0), "-1");
    }
}
