# geo-prox-rs

A really simple and fast rust library that tells whether two GPS points are close to each other by leveraging cheap lower bound of the otherwise expensive geodesic distance evaluation.

## Basic usage

```rust
use geo_types::Point; // cargo add geo-types
use geo_prox::{A_TOL, R_TOL}; // this crate
// Two nearby points.
let p1 = Point::new(23.319941, 42.698334); // (longitue, latitude) in degrees
let p2 = Point::new(23.319920, 42.698323);
assert_eq!(
    // Check if the two points are at most 15 meters apart.
    geo_prox::isclose_opt(p1, p2, 15.0, A_TOL, R_TOL, true),
    Some(true)
);
// Two far apart points.
let p1 = Point::new(23.319941, 42.698334);
let p2 = Point::new(-25.319920, -42.698323); // negative for (west, south)
assert_eq!(
    geo_prox::isclose_opt(p1, p2, 15.0, A_TOL, R_TOL, true),
    Some(false)
);
```

## Feature flags

- `geo`: Enable to resort to exact geodesic distance when the bound is loose (default disbaled).

## Benchmark

Evaluating the bound is on average 85% faster than exact geodesic distance evaluation, and the bound is *always* tight at small scale (e.g. when radius <= 1 km).
Therefore, this crate is particularly useful to detect near-duplicate location in applications handling geographical data.
