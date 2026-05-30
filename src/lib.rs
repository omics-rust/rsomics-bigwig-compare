//! Per-bin comparison of two bigWig files as a bedGraph track — deeptools
//! `bigwigCompare`.
//!
//! The genome is tiled into `bin_size`-wide bins. For each chromosome common to
//! both files the per-base bigWig values are read, averaged per bin (deeptools
//! `getCoverageFromBigwig` → `numpy.mean` over each tile slice), then combined
//! per bin by an [`Operation`] (deeptools default `log2`). Output is bedGraph
//! with adjacent equal-value bins merged (unless `--fixedStep`).
//!
//! ## Chromosome set (deeptools `writeBedGraph`)
//!
//! Only chromosomes present in BOTH files are processed, in the first file's
//! B-tree leaf order (the order pyBigWig's `chroms()` dict yields, which
//! deeptools sorts its output by). Where a chromosome's declared length differs
//! between the files, the smaller length is used.
//!
//! ## Per-bin value (deeptools `getCoverageFromBigwig`)
//!
//! Per-base values are read for `[0, chrom_len)`. With `missing_data_as_zero`
//! (deeptools default — `not --skipNonCoveredRegions`) NaN bases become 0
//! first. Each bin's value is `numpy.mean` over its `bin_size`-base slice (the
//! final bin may be shorter). Without `missing_data_as_zero`, an all-NaN bin
//! averages to NaN.
//!
//! ## Combination (deeptools `getRatio` / `compute_ratio`)
//!
//! Per bin, with `v1 = scale[0] * cov1`, `v2 = scale[1] * cov2`. If either
//! scaled value is NaN the result is NaN. Otherwise:
//!
//! - `log2` — `log2((v1 + pc0) / (v2 + pc1))`
//! - `ratio` — `(v1 + pc0) / (v2 + pc1)`
//! - `reciprocal_ratio` — `r` if `r >= 1` else `-1 / r`, where
//!   `r = (v1 + pc0) / (v2 + pc1)`
//! - `subtract` — `v1 - v2`
//! - `add` — `v1 + v2`
//! - `mean` — `(v1 + v2) / 2`
//! - `first` — `v1`
//! - `second` — `v2`
//!
//! The pseudocount `[pc0, pc1]` (deeptools default `[1, 1]`) is added only for
//! the ratio-family operations. With `--skipZeroOverZero`, a bin whose two
//! coverage values sum to zero (before any pseudocount) is dropped. bedGraph
//! values use Python's `{:g}` format.

#![allow(clippy::cast_precision_loss)]

pub mod bins;
pub mod operation;
pub mod output;

pub use operation::Operation;

use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_bbi::BigWig;
use rsomics_common::{Result, RsomicsError};

use bins::common_chroms;
use output::write_chrom;

#[derive(Debug, Clone)]
pub struct CompareOpts {
    pub bin_size: u32,
    pub operation: Operation,
    pub scale_factors: [f64; 2],
    pub pseudocount: [f64; 2],
    pub missing_data_as_zero: bool,
    pub skip_zero_over_zero: bool,
    pub fixed_step: bool,
}

impl Default for CompareOpts {
    fn default() -> Self {
        Self {
            bin_size: 50,
            operation: Operation::Log2,
            scale_factors: [1.0, 1.0],
            pseudocount: [1.0, 1.0],
            missing_data_as_zero: true,
            skip_zero_over_zero: false,
            fixed_step: false,
        }
    }
}

/// Bin both bigWigs over common chromosomes, emit bedGraph to `output`.
/// Returns the number of bedGraph lines written.
pub fn bigwig_compare(
    bw1: &Path,
    bw2: &Path,
    output: &mut dyn Write,
    opts: &CompareOpts,
) -> Result<u64> {
    let mut a = BigWig::open(bw1)?;
    let mut b = BigWig::open(bw2)?;

    let common = common_chroms(&a, &b);

    let mut out = BufWriter::with_capacity(256 * 1024, output);
    let mut lines: u64 = 0;
    for (chrom, len) in &common {
        lines += write_chrom(&mut out, &mut a, &mut b, chrom, *len, opts)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(lines)
}
