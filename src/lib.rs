//! This crate provides a fast implementation on whether two lon/lat points are
//! close to each other, based on the observation that the chord distance is a
//! tight lower bound of the geodesic distance at small scale.
//!
//! Usage:
//!
//! ```
//! use geo_types::Point;
//! let p1 = Point::new(23.319941, 42.698334);
//! let p2 = Point::new(23.319920, 42.698323);
//! assert_eq!(geo_prox::isclose_opt(p1, p2, 15.0, None, None), Some(true));
//! ```

#[cfg(feature = "geo")]
use geo::{Distance, GeodesicMeasure};
use geo_types::{CoordFloat, Point};

/// WGS84 ellipsoid constants. See https://en.wikipedia.org/wiki/Flattening for
/// details.
mod wgs84 {
    /// Semi-maojor axis (equatorial radius) in meters.
    pub const A: f64 = 6378137.0;
    /// Earth flattening.
    const F: f64 = 1.0 / 298.257223563;
    /// Eccentricity squared.
    pub const E2: f64 = F * (2.0 - F);
    /// Semi-minor axis (polar radius) in meters.
    pub const B: f64 = A * (1.0 - F);
}

/// The ECEF coordinate.
#[derive(Debug, Clone, Copy)]
struct Ecef<T: CoordFloat> {
    x: T,
    y: T,
    z: T,
}

impl From<Point<f64>> for Ecef<f64> {
    /// Convert from (longitude, latitude) in degrees to ECEF coordinate.
    fn from(pt: Point<f64>) -> Self {
        let pt = pt.to_radians();
        // Longitude and latitude in radians.
        let (lon, lat) = pt.x_y();
        // The altitude, assumed 0.0 for now.
        let alt = 0.0;
        let (lon_s, lon_c) = lon.sin_cos();
        let (lat_s, lat_c) = lat.sin_cos();
        let n = wgs84::A / (1.0 - wgs84::E2 * lat_s * lat_s).sqrt();
        let x = (n + alt) * lat_c * lon_c;
        let y = (n + alt) * lat_c * lon_s;
        let z = (n * (1.0 - wgs84::E2) + alt) * lat_s;
        Self { x, y, z }
    }
}

impl Ecef<f64> {
    /// Compute the chord distance through the earth.
    fn chord_distance(&self, other: Ecef<f64>) -> f64 {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        let dz = self.z - other.z;
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// A conservative bound on (geodesic - chord) using cubic small-distance
/// approximation derived from the spherical case.
const fn max_geodesic_minus_chord_bound(chord_dist: f64) -> f64 {
    (chord_dist * chord_dist * chord_dist) / (24.0 * wgs84::B * wgs84::B)
}

/// Test if the error `err` (non-negative) is tiny comparing with `b`
/// (non-negative).
const fn err_is_tiny(err: f64, b: f64, a_tol: f64, r_tol: f64) -> bool {
    err <= a_tol + r_tol * b
}

/// Return a lower bound and an optional upper bound of the geodesic distance
/// between two points (longitude, latitude) in degrees. If the lower bound is
/// greater than `radius`, only the lower bound will be returned; otherwise,
/// both bounds will be returned. Set `radius` to infinity if one wants to
/// ensure the upper bound is not None.
pub fn geodesic_distance_bound(
    p1: Point<f64>,
    p2: Point<f64>,
    radius: f64,
) -> (f64, Option<f64>) {
    let ecef_p1: Ecef<_> = p1.into();
    let ecef_p2: Ecef<_> = p2.into();
    let d_c = ecef_p1.chord_distance(ecef_p2);

    if d_c > radius {
        (d_c, None)
    } else {
        let err_upper_bound = max_geodesic_minus_chord_bound(d_c);
        (d_c, Some(d_c + err_upper_bound))
    }
}

/// Test whether two points (longitude, latitude) in degrees are within
/// `radius` meters in geodesic distance. Instead of running full geodesic
/// distance, we test the proximity with the chord distance lower bound. At
/// scale of a few hundred meters of `radius`, the error is tiny. Fallback
/// to [`geo::GeodesicMeasure`] if `radius` is too large to make the chord
/// distance a tight enough bound. Whether `radius` is regarded as "too
/// large" is decided by the absolute tolerance `a_tol` and relative tolerance
/// `r_tol`. Passing None to use the default tolerances, which are sensible for
/// most use cases.
///
/// # Accepted coordinate value
///
/// - longitude: -180.0..=180.0, negative for west, positive for east.
/// - latitude: -90.0..=90.0, negative for south, positive for north.
///
/// # Notes on coordinate type
///
/// We have to stick to f64 because [`geo::GeodesicMeasure`] only implements
/// [`Distance`](geo::algorithm::line_measures::Distance) on f64.
#[cfg(feature = "geo")]
pub fn isclose(
    p1: Point<f64>,
    p2: Point<f64>,
    radius: f64,
    fallback_dist: &GeodesicMeasure,
    a_tol: Option<f64>,
    r_tol: Option<f64>,
) -> bool {
    let (lb, ub_opt) = geodesic_distance_bound(p1, p2, radius);

    match ub_opt {
        None => {
            // If the lower bound already exeeds `radius`, `p1` and `p2` are
            // definitely not close.
            false
        }
        Some(ub) => {
            // The absolute tolerance, default taken from numpy.
            let a_tol = a_tol.unwrap_or(1e-8);
            // The relative tolerance, default taken from numpy.
            let r_tol = r_tol.unwrap_or(1e-5);

            if !err_is_tiny(ub - lb, ub, a_tol, r_tol) {
                // The bound is too loose. We have to resort to the exact
                // distance measure.
                fallback_dist.distance(p1, p2) <= radius
            } else {
                // Trade off a little bit precision for speed, as the bound is
                // tight enough, assuming that `a_tol` and `r_tol` values are
                // sensible.
                true
            }
        }
    }
}

/// Test whether two points (longitude, latitude) in degrees are within
/// `radius` meters in geodesic distance. Instead of running full geodesic
/// distance, we test the proximity with the chord distance lower bound. At
/// scale of a few hundred meters of `radius`, the error is tiny. Return None
/// if `radius` is too large to make the chord distance a tight enough bound.
/// Whether `radius` is regarded as "too large" is decided by the absolute
/// tolerance `a_tol` and relative tolerance `r_tol`. Passing None to use the
/// default tolerances, which are sensible for most use cases.
///
/// # Accepted coordinate value
///
/// - longitude: -180.0..=180.0, negative for west, positive for east.
/// - latitude: -90.0..=90.0, negative for south, positive for north.
///
/// # Notes on coordinate type
///
/// We have to stick to f64 because [`geo::GeodesicMeasure`] only implements
/// [`Distance`](geo::algorithm::line_measures::Distance) on f64.
pub fn isclose_opt(
    p1: Point<f64>,
    p2: Point<f64>,
    radius: f64,
    a_tol: Option<f64>,
    r_tol: Option<f64>,
) -> Option<bool> {
    let (lb, ub_opt) = geodesic_distance_bound(p1, p2, radius);

    match ub_opt {
        None => {
            // If the lower bound already exeeds `radius`, `p1` and `p2` are
            // definitely not close.
            Some(false)
        }
        Some(ub) => {
            // The absolute tolerance, default taken from numpy.
            let a_tol = a_tol.unwrap_or(1e-8);
            // The relative tolerance, default taken from numpy.
            let r_tol = r_tol.unwrap_or(1e-5);

            if !err_is_tiny(ub - lb, ub, a_tol, r_tol) {
                // The bound is too loose. We choose to abstain.
                None
            } else {
                // Trade off a little bit precision for speed, as the bound is
                // tight enough, assuming that `a_tol` and `r_tol` values are
                // sensible.
                Some(true)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use geo_types::Point;

    #[test]
    fn test_isclose() {
        let p1 = Point::new(23.319941, 42.698334);
        let p2 = Point::new(23.319920, 42.698323);
        assert_eq!(super::isclose_opt(p1, p2, 15.0, None, None), Some(true));
    }
}
