use cardio_rs::{HrvMetrics, standard_analysis};
use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;

fn analysis() {
    let path = "tests/ecg.csv";
    let signal = "ECG_Raw";
    let rate = 1000.;
    let hrv_metrics = standard_analysis!(path, signal, rate);
    black_box(hrv_metrics); // Prevents optimization that would remove the analysis entirely
}

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("standard_analysis", |b| b.iter(analysis));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
