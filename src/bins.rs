use std::collections::HashMap;

use rsomics_bbi::BigWig;
use rsomics_common::{Result, RsomicsError};

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
