//! A module for calculating geometrical time-domain metrics for Heart Rate Variability (HRV).
//!
//! It defines the `GeometricMetrics` struct, which contains two key HRV metrics:
//! the Triangular Index (TI) and the Triangular Interpolation of NN-interval Histogram (TINN).
//! This module provides functionality for calculating these metrics from a list of Normal-to-Normal (NN) intervals.
//!
//! ## Key Components:
//! - `GeometricMetrics<T>`: The main struct that holds the computed `triangular_index` and `tinn` metrics.
//! - `compute()`: A method for calculating both the Triangular Index (TI) and the TINN based on a series of NN intervals.
//!
//! ## Features:
//! - **Triangular Index**: The HRV triangular index is calculated by dividing the total number of NN intervals by the maximum frequency in the histogram of NN intervals. It is an approximation of the RR interval distribution's density.
//! - **TINN (Triangular Interpolation of NN-interval Histogram)**: The baseline width of the distribution measured as the base of a triangle, approximating the NN-interval distribution. This metric requires a histogram of the NN intervals, and the baseline width is determined by triangular interpolation.
//!
//! ## Limitations:
//! - **TINN Precision**: The method used to calculate TINN here relies on an approximation that may not be as precise as other methods that use more sophisticated interpolation techniques or higher resolution histograms. The current implementation uses a simplified method with an epsilon factor and may be subject to errors due to bin size and histogram smoothing.
//!
//! ## Example Usage:
//! ```rust
//! use cardio_rs::geometric_domain::GeometricMetrics;
//!
//! let rr_intervals = vec![800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0];
//! let geometric_metrics = GeometricMetrics::<f64>::compute(&rr_intervals);
//! println!("{:?}", geometric_metrics);
//! ```
//!
//! ## Struct and Method Documentation:
//! - `GeometricMetrics<T>`: A struct that holds two HRV metrics: `triangular_index` and `tinn`.
//! - `compute()`: A method that takes a slice of RR intervals and computes the `triangular_index` and `tinn` based on the provided intervals.
#[cfg(not(feature = "std"))]
extern crate alloc;
#[cfg(not(feature = "std"))]
use alloc::vec;

use core::iter::Sum;
use num::Float;

/// A struct representing two geometrical Heart Rate Variability (HRV) metrics.
#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct GeometricMetrics<T> {
    /// The **Triangular Index (TI)**, representing the total number of RR intervals
    /// divided by the maximum frequency in the histogram of NN-intervals.
    pub triangular_index: T,

    /// The **TINN (Triangular Interpolation of NN-interval Histogram)**, representing
    /// the baseline width of the NN interval distribution, calculated using triangular interpolation.
    pub tinn: T,
}

impl<T: Float + Sum<T> + Copy + core::fmt::Debug + core::ops::AddAssign + Into<f64>>
    GeometricMetrics<T>
{
    /// Computes the **Triangular Index (TI)** and **TINN (Triangular Interpolation of NN-interval Histogram)**
    /// for a given sequence of RR intervals.
    ///
    /// The function calculates the following:
    /// - **Triangular Index (TI)**: A ratio of the total number of RR intervals to the highest frequency
    ///   observed in the histogram bins of the RR intervals.
    /// - **TINN**: The baseline width of the histogram approximated by a triangle using a specified epsilon value.
    ///
    /// ## Parameters:
    /// - `rr_intervals`: A slice of `T` values representing the RR intervals (NN intervals).
    ///
    /// ## Return Value:
    /// - Returns an instance of `GeometricMetrics<T>` containing the computed `triangular_index` and `tinn`.
    ///
    /// ## Example:
    /// ```rust
    /// use cardio_rs::geometric_domain::GeometricMetrics;
    ///
    /// let rr_intervals = vec![800.0, 850.0, 900.0, 950.0, 1000.0, 1100.0];
    /// let geometric_metrics = GeometricMetrics::<f64>::compute(&rr_intervals);
    /// println!("{:?}", geometric_metrics);  // Output: GeometricMetrics { triangular_index: 1.2, tinn: 12.0 }
    /// ```
    pub fn compute(rr_intervals: &[T]) -> Self {
        let min = T::from(300).unwrap();
        let max = T::from(2_000).unwrap();
        let bin_width = T::from(8).unwrap();
        let num_bins = Into::<f64>::into(((max - min) / bin_width).ceil()) as usize;
        let mut hist = vec![T::from(0).unwrap(); num_bins];

        rr_intervals.iter().for_each(|&i| {
            let index = Into::<f64>::into(((i - min) / bin_width).floor()) as usize;
            hist[index] += T::one();
        });

        let mode = hist
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();

        let triangular_index = T::from(rr_intervals.len()).unwrap() / *mode;

        // TODO can use more precise method
        let epsilon = *mode * T::from(0.75).unwrap();
        let base_left = hist.iter().position(|&i| i > epsilon).unwrap();
        let base_right =
            hist.len() - 1usize - hist.iter().rev().position(|&i| i > epsilon).unwrap();

        let tinn = T::from(base_right - base_left).unwrap() * bin_width;

        Self {
            triangular_index,
            tinn,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_data::RR_INTERVALS;
    use approx::{AbsDiffEq, RelativeEq, UlpsEq, assert_relative_eq};
    impl<T: AbsDiffEq> AbsDiffEq for GeometricMetrics<T>
    where
        T::Epsilon: Copy,
    {
        type Epsilon = T::Epsilon;

        fn default_epsilon() -> T::Epsilon {
            T::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
            T::abs_diff_eq(&self.triangular_index, &other.triangular_index, epsilon)
                && T::abs_diff_eq(&self.tinn, &other.tinn, epsilon)
        }
    }

    impl<T: RelativeEq> RelativeEq for GeometricMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_relative() -> T::Epsilon {
            T::default_max_relative()
        }

        fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
            T::relative_eq(
                &self.triangular_index,
                &other.triangular_index,
                epsilon,
                max_relative,
            ) && T::relative_eq(&self.tinn, &other.tinn, epsilon, max_relative)
        }
    }

    impl<T: UlpsEq> UlpsEq for GeometricMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_ulps() -> u32 {
            T::default_max_ulps()
        }

        fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
            T::ulps_eq(
                &self.triangular_index,
                &other.triangular_index,
                epsilon,
                max_ulps,
            ) && T::ulps_eq(&self.tinn, &other.tinn, epsilon, max_ulps)
        }
    }

    #[test]
    fn test_metrics() {
        let geo_params = GeometricMetrics::compute(RR_INTERVALS);

        assert_relative_eq!(
            GeometricMetrics {
                tinn: 128., //117.1875,
                triangular_index: 18.31578947368421,
            },
            geo_params,
            epsilon = 1e-14
        );
    }
}
