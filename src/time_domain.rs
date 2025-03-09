//! A module for calculating and storing time-domain heart rate variability (HRV) metrics.
//!
//! This module defines the `TimeMetrics` struct, which holds various HRV parameters derived from RR intervals.
//! These time-domain metrics are essential for evaluating the autonomic nervous system's regulation of the heart and understanding heart rate variability over time.
//!
//! The key metrics include:
//! - `RMSSD`: Root Mean Square of Successive Differences
//! - `SDNN`: Standard Deviation of NN Intervals
//! - `PNN50`: Percentage of Successive RR Intervals > 50 ms
//! - `Mean HR`: Mean Heart Rate
//! - `SDSD`: Standard Deviation of Successive Differences
//! - `AVNN`: Average NN Interval
//! - `CVSD`: Coefficient of Variation of Successive Differences
//!
//! The `TimeMetrics` struct allows easy access to these metrics after calculating them from a set of RR intervals.
//! This is useful for cardiovascular health monitoring, stress assessment, and various other applications in physiology and healthcare analytics.
//!
//! # `no_std` Compatibility
//!
//! This module is designed to work in both standard (`std`) and `no_std` environments. It does not rely on the Rust standard library,
//! making it suitable for use in embedded systems or other environments where `std` is unavailable.
//!
//! # Example
//! ```
//! use time_metrics::TimeMetrics;
//!
//! let rr_intervals = vec![800, 810, 780, 850, 900];
//! let metrics = TimeMetrics::compute(&rr_intervals);
//!
//! println!("{:?}", metrics);
//! ```
pub mod time_metrics {
    // The TimeMetrics struct and the compute function would be inside this module.
    // The rest of the code goes here...
}

extern crate alloc;
use alloc::vec::Vec;
use core::iter::Sum;
use num::Float;

/// A struct representing various heart rate variability (HRV) time-domain metrics.
///
/// The `TimeMetrics` struct holds key HRV parameters computed from RR intervals.
/// These metrics are essential for assessing the short-term and long-term variability of heart rate, which can provide insights into the autonomic nervous system's regulation of the heart.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct TimeMetrics<T> {
    /// Root Mean Square of Successive Differences (RMSSD).
    /// A measure of short-term HRV derived from the differences between successive RR intervals.
    pub rmssd: T,

    /// Standard Deviation of NN Intervals (SDNN).
    /// A measure of overall HRV that reflects the variability in RR intervals over time.
    pub sdnn: T,

    /// Percentage of Successive RR Intervals > 50 ms (PNN50).
    /// This metric quantifies the proportion of intervals that differ by more than 50 ms, indicating parasympathetic nervous system activity.
    pub pnn50: T,

    /// Mean Heart Rate (mean_hr).
    /// The average heart rate calculated from the RR intervals over the recording period.
    pub mean_hr: T,

    /// Standard Deviation of Successive Differences (SDSD).
    /// A short-term HRV metric that measures the standard deviation of the differences between successive RR intervals.
    pub sdsd: T,

    /// Average NN Interval (AVNN).
    /// The mean of the RR intervals, reflecting the average time between heartbeats.
    pub avnn: T,

    /// Coefficient of Variation of Successive Differences (CVSD).
    /// A normalized measure of short-term HRV, calculated as the ratio of the standard deviation to the mean of the successive RR interval differences.
    pub cvsd: T,
}

impl<T: Float + Sum<T> + Copy + core::fmt::Debug> TimeMetrics<T> {
    /// Computes various time-domain HRV metrics from the given RR intervals.
    ///
    /// The function calculates the following HRV metrics:
    /// - **SDNN (Standard Deviation of Normal-to-Normal intervals)**: A measure of overall HRV, indicating long-term variability.
    /// - **RMSSD (Root Mean Square of Successive Differences)**: A short-term HRV measure that reflects the variability between successive RR intervals.
    /// - **SDSD (Standard Deviation of Successive Differences)**: Another short-term HRV measure based on the standard deviation of successive RR interval differences.
    /// - **PNN50 (Percentage of Successive RR Interval Differences > 50 ms)**: A measure of parasympathetic nervous system activity.
    /// - **Mean HR (Mean Heart Rate)**: The average heart rate, derived from RR intervals.
    /// - **AVNN (Average NN Interval)**: The mean of all RR intervals.
    /// - **CVSD (Coefficient of Variation of Successive Differences)**: The ratio of RMSSD to the average RR interval, indicating relative variability.
    ///
    /// # Arguments
    ///
    /// * `rr_intervals` - A slice of `RR` intervals representing the time between successive heartbeats (in milliseconds).
    ///
    /// # Returns
    ///
    /// Returns a `TimeMetrics` struct containing the computed HRV metrics.
    ///
    /// # Example
    ///
    /// ```
    /// let rr_intervals = vec![800.0, 820.0, 810.0, 780.0, 790.0];
    /// let metrics = TimeMetrics::compute(&rr_intervals);
    /// println!("{:?}", metrics);
    /// ```
    pub fn compute(rr_intervals: &[T]) -> Self {
        let rr_diffs: Vec<T> = rr_intervals.windows(2).map(|i| i[1] - i[0]).collect();
        let rr_intervals_mean: T =
            rr_intervals.iter().copied().sum::<T>() / T::from(rr_intervals.len()).unwrap();
        let rr_diffs_mean: T =
            rr_diffs.iter().copied().sum::<T>() / T::from(rr_diffs.len()).unwrap();

        let variance = rr_intervals
            .iter()
            .map(|&x| (x - rr_intervals_mean) * (x - rr_intervals_mean))
            .sum::<T>()
            / T::from(rr_diffs.len()).unwrap();

        let sdnn = variance.sqrt();

        let rmssd =
            (rr_diffs.iter().map(|&x| x * x).sum::<T>() / T::from(rr_diffs.len()).unwrap()).sqrt();

        let sdsd = (rr_diffs
            .iter()
            .map(|&x| (x - rr_diffs_mean) * (x - rr_diffs_mean))
            .sum::<T>()
            / T::from(rr_diffs.len()).unwrap())
        .sqrt();

        let pnn50 = T::from(
            rr_diffs
                .iter()
                .filter(|&&x| x.abs() > T::from(50.0).unwrap())
                .count(),
        )
        .unwrap()
            / T::from(rr_diffs.len()).unwrap();

        let mean_hr = rr_intervals
            .iter()
            .map(|&i| T::from(60_000).unwrap() / i)
            .sum::<T>()
            / T::from(rr_intervals.len()).unwrap();

        Self {
            sdnn,
            rmssd,
            pnn50,
            mean_hr,
            sdsd,
            avnn: rr_intervals_mean,
            cvsd: rmssd / rr_intervals_mean,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::test_data::RR_INTERVALS;
    use approx::AbsDiffEq;
    use approx::RelativeEq;
    use approx::UlpsEq;
    use approx::assert_relative_eq;
    impl<T: AbsDiffEq> AbsDiffEq for TimeMetrics<T>
    where
        T::Epsilon: Copy,
    {
        type Epsilon = T::Epsilon;

        fn default_epsilon() -> T::Epsilon {
            T::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
            T::abs_diff_eq(&self.rmssd, &other.rmssd, epsilon)
                && T::abs_diff_eq(&self.sdnn, &other.sdnn, epsilon)
                && T::abs_diff_eq(&self.pnn50, &other.pnn50, epsilon)
                && T::abs_diff_eq(&self.mean_hr, &other.mean_hr, epsilon)
                && T::abs_diff_eq(&self.sdsd, &other.sdsd, epsilon)
                && T::abs_diff_eq(&self.cvsd, &other.cvsd, epsilon)
        }
    }

    impl<T: RelativeEq> RelativeEq for TimeMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_relative() -> T::Epsilon {
            T::default_max_relative()
        }

        fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
            T::relative_eq(&self.rmssd, &other.rmssd, epsilon, max_relative)
                && T::relative_eq(&self.sdnn, &other.sdnn, epsilon, max_relative)
                && T::relative_eq(&self.pnn50, &other.pnn50, epsilon, max_relative)
                && T::relative_eq(&self.mean_hr, &other.mean_hr, epsilon, max_relative)
                && T::relative_eq(&self.sdsd, &other.sdsd, epsilon, max_relative)
                && T::relative_eq(&self.cvsd, &other.cvsd, epsilon, max_relative)
        }
    }

    impl<T: UlpsEq> UlpsEq for TimeMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_ulps() -> u32 {
            T::default_max_ulps()
        }

        fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
            T::ulps_eq(&self.rmssd, &other.rmssd, epsilon, max_ulps)
                && T::ulps_eq(&self.sdnn, &other.sdnn, epsilon, max_ulps)
                && T::ulps_eq(&self.pnn50, &other.pnn50, epsilon, max_ulps)
                && T::ulps_eq(&self.mean_hr, &other.mean_hr, epsilon, max_ulps)
                && T::ulps_eq(&self.sdsd, &other.sdsd, epsilon, max_ulps)
                && T::ulps_eq(&self.cvsd, &other.cvsd, epsilon, max_ulps)
        }
    }

    #[test]
    fn test_metrics() {
        let time_params = TimeMetrics::compute(&RR_INTERVALS);

        assert_relative_eq!(
            TimeMetrics {
                rmssd: 59.99807873965272,
                sdnn: 69.54843274772418,
                pnn50: 0.414985590778098,
                mean_hr: 70.46146532393496,
                sdsd: 59.99789159342577,
                cvsd: 0.06999675282912315,
                avnn: 857.1551724137931,
            },
            time_params,
            epsilon = 1e-14
        );
    }
}
