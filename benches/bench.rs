use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bigwig_compare(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bigwig-compare");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let a = manifest.join("tests/golden/a.bw");
    let b = manifest.join("tests/golden/b.bw");
    let out = tempfile::NamedTempFile::new().unwrap();

    c.bench_function("rsomics-bigwig-compare golden", |bench| {
        bench.iter(|| {
            let status = Command::new(black_box(bin))
                .args([
                    "--bigwig1",
                    a.to_str().unwrap(),
                    "--bigwig2",
                    b.to_str().unwrap(),
                    "--operation",
                    "log2",
                    "-o",
                    out.path().to_str().unwrap(),
                ])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_bigwig_compare);
criterion_main!(benches);
