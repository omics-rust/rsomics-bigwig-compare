//! Chromosome intersection and per-bin averaging from bigWig.

use std::collections::HashMap;

use rsomics_bbi::BigWig;
use rsomics_common::{Result, RsomicsError};

/// Chromosomes present in both files, in the first file's order.
/// Length is the smaller of the two when they disagree (deeptools' `min`).
pub(crate) fn common_chroms(a: &BigWig, b: &BigWig) -> Vec<(String, u32)> {
    let b_lens: HashMap<&str, u32> = b.chroms().collect();
    a.chroms()
        .filter_map(|(name, alen)| {
            b_lens
                .get(name)
                .map(|&blen| (name.to_string(), alen.min(blen)))
        })
        .collect()
}

/// Per-bin average over `[0, len)`, mirroring deeptools `getCoverageFromBigwig`.
///
/// Per-base values are NaN→0 when `missing_as_zero`; each bin is a
/// `bin_size`-base tile (last tile may be shorter), averaged with numpy.mean
/// semantics (an all-NaN tile without `missing_as_zero` yields NaN).
pub(crate) fn binned_values(
    bw: &mut BigWig,
    chrom: &str,
    len: u32,
    bin_size: u32,
    missing_as_zero: bool,
) -> Result<Vec<f64>> {
    let per_base = bw.values(chrom, 0, len)?.ok_or_else(|| {
        RsomicsError::InvalidInput(format!("chromosome {chrom} vanished from bigWig"))
    })?;

    let n_bins = (len as usize).div_ceil(bin_size as usize);
    let bs = bin_size as usize;
    let mut bins = Vec::with_capacity(n_bins);
    let mut x = 0usize;
    while x < per_base.len() {
        let end = (x + bs).min(per_base.len());
        let slice = &per_base[x..end];
        // deeptools zero-fills NaN bases before the mean when missing_as_zero;
        // otherwise numpy.mean of a slice containing any NaN is itself NaN.
        let mut sum = 0.0f64;
        let mut any_nan = false;
        for &v in slice {
            if v.is_nan() {
                any_nan = true;
            } else {
                sum += f64::from(v);
            }
        }
        let mean = if slice.is_empty() || (any_nan && !missing_as_zero) {
            f64::NAN
        } else {
            sum / slice.len() as f64
        };
        bins.push(mean);
        x += bs;
    }
    Ok(bins)
}
