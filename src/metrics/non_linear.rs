//! A module for computing non-linear heart rate variability (HRV) metrics.
//!
//! This module defines the `NonLinearMetrics` struct, which holds various non-linear HRV parameters derived from RR intervals.
//! These non-linear metrics are essential for assessing the complexity and fractal properties of heart rate variability, providing insights into the autonomic nervous system's regulation of the heart.
//!
//! The key metrics include:
//! - `SampEn`: Sample Entropy
//! - `DFA`: Detrended Fluctuation Analysis
//! - `LZC`: Lempel-Ziv Complexity
//!
//! The `NonLinearMetrics` struct allows easy access to these metrics after calculating them from a set of RR intervals.
//! These metrics are valuable for understanding the non-linear dynamics of heart rate variability, which can be indicative of cardiovascular health, stress levels, and overall autonomic nervous system function.
//!
//! # Example Usage
//! ```
//! use cardio_rs::non_linear::NonLinearMetrics;
//!
//! let rr_intervals = vec![800.0, 820.0, 810.0, 780.0, 790.0];
//! let metrics = NonLinearMetrics::compute_default(&rr_intervals);
//!
//! println!("{:?}", metrics);
//! ```

#[cfg(not(feature = "std"))]
extern crate alloc;
#[cfg(not(feature = "std"))]
use alloc::{vec, vec::Vec};
use core::iter::Sum;
use num::Float;

/// A struct representing various non-linear heart rate variability (HRV) metrics.
///
/// The `NonLinearMetrics` struct holds key HRV parameters computed from RR intervals using non-linear techniques.
/// These metrics provide insight into the complexity and fractal characteristics of heart rate dynamics, which reflect autonomic nervous system regulation of the heart.
/// These non-linear metrics are important for understanding the self-organizing behavior of the cardiovascular system over time.
#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct NonLinearMetrics<T> {
    /// Sample Entropy (SampEn).
    /// A measure of the regularity and complexity of the RR intervals. Lower values indicate a more regular signal.
    pub samp_en: T,

    /// Detrended Fluctuation Analysis (DFA).
    /// A measure of long-range correlations in the RR intervals, indicating fractal-like properties of the heart rate variability signal.
    pub dfa: T,

    /// Lempel-Ziv Complexity (LZC).
    /// A measure of the signal's complexity based on the number of distinct patterns present in the RR intervals.
    pub lzc: T,
}

impl<T: Float + Sum<T> + Copy + core::fmt::Debug + core::ops::AddAssign> NonLinearMetrics<T> {
    /// Computes various non-linear HRV metrics from the given RR intervals with default parameters.
    ///
    /// The function calculates the following HRV metrics:
    /// - **SampEn (Sample Entropy)**: A measure of regularity and complexity of the signal.
    /// - **DFA (Detrended Fluctuation Analysis)**: A measure of long-range correlations in the RR intervals.
    /// - **LZC (Lempel-Ziv Complexity)**: A measure of the signal's complexity based on the number of distinct patterns.
    ///
    /// # Arguments
    ///
    /// * `rr_intervals` - A slice of RR intervals representing the time between successive heartbeats (in milliseconds).
    ///
    /// # Returns
    ///
    /// Returns a `NonLinearMetrics` struct containing the computed non-linear HRV metrics.
    pub fn compute_default(rr_intervals: &[T]) -> Self {
        Self {
            samp_en: Self::compute_samp_en(rr_intervals, &T::from(0.2).unwrap()),
            dfa: Self::compute_dfa(rr_intervals, None),
            lzc: Self::compute_lzc(rr_intervals),
        }
    }

    /// Computes various non-linear HRV metrics from the given RR intervals with custom tolerance and scales.
    ///
    /// The function calculates the following HRV metrics:
    /// - **SampEn (Sample Entropy)**: A measure of regularity and complexity of the signal.
    /// - **DFA (Detrended Fluctuation Analysis)**: A measure of long-range correlations in the RR intervals.
    /// - **LZC (Lempel-Ziv Complexity)**: A measure of the signal's complexity based on the number of distinct patterns.
    ///
    /// # Arguments
    ///
    /// * `rr_intervals` - A slice of RR intervals representing the time between successive heartbeats (in milliseconds).
    /// * `tolerance` - A tolerance value used for calculating sample entropy.
    /// * `scales` - An optional list of scale sizes to be used for DFA calculation.
    ///
    /// # Returns
    ///
    /// Returns a `NonLinearMetrics` struct containing the computed non-linear HRV metrics.
    pub fn compute(rr_intervals: &[T], tolerance: &T, scales: Option<&[usize]>) -> Self {
        Self {
            samp_en: Self::compute_samp_en(rr_intervals, tolerance),
            dfa: Self::compute_dfa(rr_intervals, scales),
            lzc: Self::compute_lzc(rr_intervals),
        }
    }

    /// Computes the Sample Entropy (SampEn) from the given RR intervals and tolerance value.
    ///
    /// SampEn is a measure of the complexity and regularity of the time series.
    /// Lower values indicate a more regular signal, while higher values suggest more complexity.
    ///
    /// # Arguments
    ///
    /// * `rr_intervals` - A slice of RR intervals representing the time between successive heartbeats.
    /// * `tolerance` - A tolerance value used to define how similar two sequences must be to be considered similar.
    ///
    /// # Returns
    ///
    /// Returns the Sample Entropy (SampEn) of the RR intervals as a floating point value.
    pub fn compute_samp_en(rr_intervals: &[T], tolerance: &T) -> T {
        let rr_intervals_mean: T =
            rr_intervals.iter().copied().sum::<T>() / T::from(rr_intervals.len()).unwrap();

        let variance = rr_intervals
            .iter()
            .map(|&x| (x - rr_intervals_mean) * (x - rr_intervals_mean))
            .sum::<T>()
            / T::from(rr_intervals.len() - 1).unwrap();

        let r = *tolerance * variance.sqrt();

        let pairs: Vec<(T, T, T)> = rr_intervals
            .windows(3)
            .map(|i| (i[0], i[1], i[2]))
            .collect();

        let mut count_a = 0usize;
        let mut count_b = 0usize;
        for (i, j) in pairs.iter().enumerate() {
            for l in pairs.iter().skip(i + 1) {
                if i == 0 && T::max((j.0 - l.0).abs(), (j.1 - l.1).abs()) < r {
                    count_b += 1;
                }

                if T::max((j.2 - l.2).abs(), (j.1 - l.1).abs()) < r {
                    count_b += 1;
                }
                if T::max(
                    T::max((j.0 - l.0).abs(), (j.1 - l.1).abs()),
                    (j.2 - l.2).abs(),
                ) < r
                {
                    count_a += 1;
                }
            }
        }
        let b = T::from(count_b).unwrap() / T::from(rr_intervals.len() - 1).unwrap();
        let a = T::from(count_a).unwrap() / T::from(rr_intervals.len() - 2).unwrap();

        (b / a).ln()
    }

    /// Performs linear regression on the given x and y values and returns the slope and intercept.
    ///
    /// # Arguments
    ///
    /// * `x` - A slice of x-values.
    /// * `y` - A slice of y-values.
    ///
    /// # Returns
    ///
    /// Returns a tuple `(a, b)`, where `a` is the slope and `b` is the intercept of the regression line.
    fn linear_regression(x: &[T], y: &[T]) -> (T, T) {
        let n = T::from(y.len()).unwrap();
        let x_sum = x.iter().copied().sum::<T>();
        let y_sum = y.iter().copied().sum::<T>();
        let xy_sum = y.iter().zip(x.iter()).map(|(&y, &x)| x * y).sum::<T>();
        let xx_sum = x.iter().map(|&i| i * i).sum::<T>();

        let a = (n * xy_sum - x_sum * y_sum) / (n * xx_sum - x_sum * x_sum);
        let b = (y_sum - a * x_sum) / n;

        (a, b)
    }

    /// Computes the Detrended Fluctuation Analysis (DFA) from the given RR intervals and scale values.
    ///
    /// DFA measures the long-range correlations in the RR intervals, indicating fractal-like properties of heart rate variability.
    /// This method applies linear regression to segments of the data and calculates the fluctuation at each scale.
    ///
    /// # Arguments
    ///
    /// * `rr_intervals` - A slice of RR intervals representing the time between successive heartbeats.
    /// * `scales` - An optional slice of scale sizes to use for DFA. If None, scales are generated automatically.
    ///
    /// # Returns
    ///
    /// Returns the scaling exponent (α) from DFA, indicating the long-term correlation.
    pub fn compute_dfa(rr_intervals: &[T], scales: Option<&[usize]>) -> T {
        let rr_intervals_mean: T =
            rr_intervals.iter().copied().sum::<T>() / T::from(rr_intervals.len()).unwrap();

        let y: Vec<T> = rr_intervals
            .iter()
            .map(|&i| (i - rr_intervals_mean))
            .scan(T::from(0).unwrap(), |s, i| {
                *s += i;
                Some(*s)
            })
            .collect();

        let scales = if let Some(scales) = scales {
            scales
        } else {
            let min_win_size = 4;
            let max_win_size = rr_intervals.len() / 10;
            &core::iter::successors(Some(min_win_size), |&prev| {
                let next = (prev as f64 * 1.25).round() as usize;
                if next < max_win_size {
                    Some(next)
                } else {
                    None
                }
            })
            .collect::<Vec<usize>>()
        };

        let mut fluctuations: Vec<T> = vec![];

        for win_size in scales {
            let segments = y.chunks_exact(*win_size);
            let x: Vec<T> = (0..*win_size).map(|i| T::from(i).unwrap()).collect();

            let f = segments.map(|i| {
                let (a, b) = Self::linear_regression(&x, i);
                let i = i
                    .iter()
                    .enumerate()
                    .map(move |(j, &l)| l - a * T::from(j).unwrap() - b);
                let n = T::from(i.len()).unwrap();
                let mean = i.clone().sum::<T>() / n;
                let i = i.map(|i| (i - mean) * (i - mean));
                (i.sum::<T>() / n).sqrt()
            });
            fluctuations.push(f.clone().sum::<T>() / T::from(f.len()).unwrap());
        }

        let (a, _) = Self::linear_regression(
            &scales
                .iter()
                .map(|i| T::from(*i).unwrap().log(T::from(10).unwrap()))
                .collect::<Vec<T>>(),
            &fluctuations
                .iter()
                .map(|i| i.log(T::from(10).unwrap()))
                .collect::<Vec<T>>(),
        );
        a
    }

    /// Computes the complexity of a binary sequence (naive implementation).
    ///
    /// This is used as part of the Lempel-Ziv Complexity (LZC) calculation. It counts the distinct patterns in a binary sequence.
    ///
    /// # Arguments
    ///
    /// * `binary_sequence` - A slice of binary values (0 or 1).
    ///
    /// # Returns
    ///
    /// Returns the complexity of the binary sequence as a floating point value.
    fn naive_complexity(binary_sequence: &[u8]) -> T {
        let length = binary_sequence.len();
        let mut complexity = 1;
        let mut u = 0;
        let mut v = 1;
        let mut w = 1;
        let mut v_max = 1;

        while w <= length {
            if binary_sequence[u + v - 1..u + v] == binary_sequence[w + v - 1..w + v] {
                v += 1;
                if w + v >= length {
                    complexity += 1;
                    break;
                }
            } else {
                if v > v_max {
                    v_max = v;
                }
                u += 1;

                if u == w {
                    complexity += 1;

                    w += v_max;
                    if w > length {
                        break;
                    } else {
                        u = 0;
                        v = 1;
                        v_max = 1;
                    }
                } else {
                    v = 1;
                }
            }
        }
        T::from(complexity).unwrap()
    }

    /// Computes the Lempel-Ziv Complexity (LZC) of the RR intervals.
    ///
    /// LZC is a measure of the complexity of a binary sequence generated from the RR intervals.
    /// It quantifies the number of distinct patterns that appear in the binary sequence.
    ///
    /// # Arguments
    ///
    /// * `rr_intervals` - A slice of RR intervals representing the time between successive heartbeats.
    ///
    /// # Returns
    ///
    /// Returns the Lempel-Ziv Complexity (LZC) as a floating point value.
    pub fn compute_lzc(rr_intervals: &[T]) -> T {
        let rr_intervals_mean: T =
            rr_intervals.iter().copied().sum::<T>() / T::from(rr_intervals.len()).unwrap();
        let sequences = rr_intervals
            .iter()
            .map(|i| if *i <= rr_intervals_mean { 0 } else { 1 })
            .collect::<Vec<u8>>();
        let complexity = Self::naive_complexity(&sequences);
        (complexity * T::from(sequences.len()).unwrap().log(T::from(2).unwrap()))
            / T::from(sequences.len()).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_data::RR_INTERVALS;
    use approx::assert_relative_eq;

    #[test]
    fn test_samp_en_default() {
        let samp_en = NonLinearMetrics::compute_samp_en(RR_INTERVALS, &0.2);
        assert_relative_eq!(samp_en, 1.779461941417649, max_relative = 0.005,);
    }

    #[test]
    fn test_samp_en() {
        let samp_en = NonLinearMetrics::compute_samp_en(RR_INTERVALS, &0.5);
        assert_relative_eq!(samp_en, 0.9064624213913623, max_relative = 0.005,);
    }

    #[test]
    fn test_regression() {
        let a = NonLinearMetrics::linear_regression(RR_INTERVALS, RR_INTERVALS);
        assert_eq!(a, (1., 0.));

        let a = NonLinearMetrics::linear_regression(&[0., 1., 2., 3., 4.], &[1., 2., 3., 4., 5.]);
        assert_eq!(a, (1., 1.));
    }

    #[test]
    fn test_dfa_default() {
        let dfa = NonLinearMetrics::compute_dfa(RR_INTERVALS, None);
        // Reference is a log scale, we have approximative log scale
        assert_relative_eq!(dfa, 1.1133025986116414, max_relative = 0.2,);
    }

    #[test]
    fn test_dfa_custom() {
        let dfa = NonLinearMetrics::compute_dfa(RR_INTERVALS, Some(&[4, 8, 16]));
        assert_relative_eq!(dfa, 1.362119197141086, max_relative = 0.05,);

        let dfa = NonLinearMetrics::compute_dfa(RR_INTERVALS, Some(&[4, 8, 16, 64]));
        assert_relative_eq!(dfa, 0.9100673766882059, max_relative = 0.05,);
    }

    #[test]
    fn test_complexity() {
        let sequences = &[1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0];
        let c = NonLinearMetrics::<f64>::naive_complexity(sequences);
        assert_eq!(c, 6.);
    }

    #[test]
    fn test_lzc() {
        let lzc = NonLinearMetrics::compute_lzc(RR_INTERVALS);
        assert_relative_eq!(lzc, 0.679317292769438, epsilon = 1e-6);
    }
}
