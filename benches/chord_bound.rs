//! Benchmarkd the speed of the chord distance bound comparing with the exact
//! geodesic distance.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use geo::{Destination, Distance, Geodesic, Point};

fn distance_baseline(p: Point<f64>, q: Point<f64>) -> f64 {
    Geodesic.distance(p, q)
}

fn chord_bound(p: Point<f64>, q: Point<f64>) -> (f64, f64) {
    geo_prox::geodesic_distance_bound(p, q)
}

fn linspace(min: f64, max: f64, num: u32) -> Vec<f64> {
    let step = (max - min) / (num - 1) as f64;
    (0..num).map(move |i| min + (i as f64) * step).collect()
}

fn bench_chord_bound(c: &mut Criterion) {
    let mut group = c.benchmark_group("chord_bound");
    for log_distance in linspace(-2.0, 7.0, 10) {
        let distance = 10f64.powf(log_distance);
        // An example bearing.
        let bearing = 90f64;
        // An example source point.
        let p = Point::new(45f64, 45f64);
        let q = Geodesic.destination(p, bearing, distance);
        group.bench_with_input(
            BenchmarkId::new("exact", log_distance),
            &(p, q),
            |b, (p, q)| b.iter(|| distance_baseline(*p, *q)),
        );
        group.bench_with_input(
            BenchmarkId::new("chord_bound", log_distance),
            &(p, q),
            |b, (p, q)| b.iter(|| chord_bound(*p, *q)),
        );
    }
    group.finish();
}

criterion_group!(benches, bench_chord_bound);
criterion_main!(benches);
