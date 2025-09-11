//! Run bound validity test on grid over the globe.

use geo::{Destination, Geodesic, Point};

fn linspace(min: f64, max: f64, num: u32) -> Vec<f64> {
    let step = (max - min) / (num - 1) as f64;
    (0..num).map(move |i| min + (i as f64) * step).collect()
}

#[test]
fn grid_test() {
    let tol = 1e-8; // Floating point error tolerance.
    let lons = linspace(-180.0, 180.0, 37);
    let lats = linspace(-90.0, 90.0, 19);
    for lon in lons.iter().copied() {
        for lat in lats.iter().copied() {
            let src = Point::new(lon, lat);
            for bearing in linspace(0.0, 360.0, 36) {
                for log_distance in linspace(-2.0, 7.0, 100) {
                    let distance = 10f64.powf(log_distance);
                    let dst = Geodesic.destination(src, bearing, distance);
                    let (lb, ub) = geo_prox::geodesic_distance_bound(src, dst);
                    assert!(
                        0.0 <= lb
                            && lb <= distance + tol
                            && distance <= ub + tol,
                        "lb: {} | distance: {} | ub: {} @ ({}, {}) ~ ({}, {})",
                        lb,
                        distance,
                        ub,
                        src.x(),
                        src.y(),
                        dst.x(),
                        dst.y()
                    );
                }
            }
        }
    }
}
