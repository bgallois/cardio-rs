//! A module for computing heart rate variability (HRV) frequency-domain metrics.
//!
//! This module provides functionality to compute frequency-domain HRV metrics from a series of RR intervals, which are the time intervals between successive heartbeats.
//! These metrics—Low Frequency (LF), High Frequency (HF), and Very Low Frequency (VLF)—are commonly used to analyze autonomic nervous system function and assess heart health.
//!
//! The module includes two primary functions:
//! - `compute_sampled`: Computes HRV frequency-domain metrics based on a series of pre-sampled RR intervals.
//! - `compute`: Computes HRV frequency-domain metrics by first interpolating raw RR intervals to a uniform sampling rate and then calculating spectral density.
//!
//! The HRV frequency-domain metrics provide insight into the balance between the sympathetic and parasympathetic branches of the autonomic nervous system. LF and HF are often used as markers of parasympathetic and sympathetic activity, respectively. VLF is associated with long-term variability and regulatory processes.
//!
//! # Example Usage
//! ```
//! use cardio_rs::frequency_domain::FrequencyMetrics;
//! use cardio_rs::utils::test_data::RR_INTERVALS;
//!
//! let rate = 4.0;
//! let frequency_metrics = FrequencyMetrics::compute_sampled(RR_INTERVALS, rate);
//! ```
//!
//! # Notes
//!
//! - LF: Typically ranges from 0.04 Hz to 0.15 Hz and reflects the balance between the sympathetic and parasympathetic nervous systems.
//! - HF: Typically ranges from 0.15 Hz to 0.40 Hz and is often associated with parasympathetic nervous system activity (specifically respiratory sinus arrhythmia).
//! - VLF: Typically ranges from 0.003 Hz to 0.04 Hz and reflects long-term regulatory processes.
//!
//! The module uses Welch's periodogram method for estimating the power spectral density (PSD) and applies trapezoidal integration to compute the total power within each frequency band.
//!
//! **Note:** The computed values heavily depend on the specific computational techniques used.
//! This implementation provides results that are within a few percent of the Python `hrv-analysis` package but may differ slightly from those obtained using `NeuroKit2` with a 4Hz resampling rate.

#[cfg(not(feature = "std"))]
extern crate alloc;
use crate::utils::welch::{Periodogram, WelchBuilder};
#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
use num::Float;

/// A struct representing frequency-domain heart rate variability (HRV) metrics.
///
/// The `FrequencyMetrics` struct holds key frequency-domain parameters that are calculated from heart rate variability (HRV) data.
/// These metrics are derived by analyzing the power spectral density (PSD) of the RR intervals, which can provide insights into the autonomic nervous system's regulation of the heart and overall cardiovascular health.
///
/// These frequency-domain metrics provide a deeper understanding of heart rate variability by focusing on the different frequency bands that reflect various physiological processes and autonomic nervous system regulation.
#[derive(Debug, PartialEq, Clone, Copy, Default)]
pub struct FrequencyMetrics<T> {
    /// Low Frequency (LF) power. Reflects sympathetic and parasympathetic balance, typically between 0.04 and 0.15 Hz.
    pub lf: T,

    /// High Frequency (HF) power. Primarily associated with parasympathetic nervous system activity, typically between 0.15 and 0.40 Hz.
    pub hf: T,

    /// Very Low Frequency (VLF) power. Represents the power of frequencies below 0.04 Hz, associated with long-term regulatory processes.
    pub vlf: T,
}

impl<
    T: Float
        + core::iter::Sum<T>
        + Copy
        + core::fmt::Debug
        + num::Signed
        + 'static
        + core::ops::AddAssign
        + core::marker::Send
        + core::marker::Sync
        + num::FromPrimitive,
> FrequencyMetrics<T>
{
    fn trapezoidal(iter: impl Iterator<Item = (T, T)>) -> T {
        let mut sum = T::zero();
        let mut prev_x = None;
        let mut prev_y = None;

        for (x, y) in iter {
            if let Some((prev_x_value, prev_y_value)) = prev_x.zip(prev_y) {
                let dx = x - prev_x_value;
                sum += T::from(0.5).unwrap() * (prev_y_value + y) * dx;
            }
            prev_x = Some(x);
            prev_y = Some(y);
        }

        sum
    }

    fn interpolate(x: &[T], y: &[T], xp: impl Iterator<Item = T>) -> Vec<T> {
        let m: Vec<T> = x
            .windows(2)
            .zip(y.windows(2))
            .map(|(i, j)| (j[1] - j[0]) / (i[1] - i[0]))
            .collect();
        xp.map(|i| {
            let xi = x.partition_point(|&p| p < i).saturating_sub(1);
            y[xi] + m[xi] * (i - x[xi])
        })
        .collect()
    }

    /// Computes frequency-domain metrics (LF, HF, VLF) from sampled RR intervals.
    ///
    /// The `compute_sampled` function computes the frequency-domain metrics of heart rate variability (HRV) using the Welch method for spectral density estimation.
    /// It requires a list of sampled RR intervals and a specified sampling rate to estimate the power spectral density (PSD) and extract the frequency components in the low-frequency (LF), high-frequency (HF), and very low-frequency (VLF) bands.
    ///
    /// This method employs Welch's periodogram to calculate the spectral density, with specified frequency bands for LF, HF, and VLF, and computes the total power in these bands using trapezoidal integration.
    ///
    /// # Parameters
    /// - `sampled_rr_intervals`: A slice of sampled RR intervals, representing the time differences between successive heartbeats, after being adjusted for the mean interval.
    /// - `rate`: The sampling rate used to interpolate the RR intervals, expressed in Hz.
    ///
    /// # Returns
    /// A `FrequencyMetrics` struct containing the computed LF, HF, and VLF values, which represent the power in their respective frequency bands.
    ///
    /// # Example
    /// ```
    /// use cardio_rs::frequency_domain::FrequencyMetrics;
    /// use cardio_rs::utils::test_data::RR_INTERVALS;
    ///
    /// let rate = 4.0;
    /// let frequency_metrics = FrequencyMetrics::compute_sampled(RR_INTERVALS, rate);
    /// ```
    ///
    /// # Notes
    /// The frequency bands for LF, HF, and VLF are typically:
    /// - LF: 0.04 - 0.15 Hz
    /// - HF: 0.15 - 0.40 Hz
    /// - VLF: 0.003 - 0.04 Hz
    pub fn compute_sampled(sampled_rr_intervals: &[T], rate: T) -> Self {
        let mean = sampled_rr_intervals.iter().copied().sum::<T>()
            / T::from_usize(sampled_rr_intervals.len()).unwrap();
        let sampled_rr_intervals: Vec<T> = sampled_rr_intervals.iter().map(|&i| i - mean).collect();
        let welch = WelchBuilder::new(sampled_rr_intervals)
            .with_fs(rate)
            .with_overlap_size(128)
            .with_segment_size(256)
            .build();

        let lf = welch
            .frequencies()
            .zip(welch.periodogram())
            .filter_map(|(f, psd_value)| {
                if f >= T::from(0.004).unwrap() && f < T::from(0.15).unwrap() {
                    Some((f, psd_value))
                } else {
                    None
                }
            });

        let hf = welch
            .frequencies()
            .zip(welch.periodogram())
            .filter_map(|(f, psd_value)| {
                if f >= T::from(0.15).unwrap() && f < T::from(0.4).unwrap() {
                    Some((f, psd_value))
                } else {
                    None
                }
            });

        let vlf = welch
            .frequencies()
            .zip(welch.periodogram())
            .filter_map(|(f, psd_value)| {
                if f >= T::from(0.003).unwrap() && f < T::from(0.04).unwrap() {
                    Some((f, psd_value))
                } else {
                    None
                }
            });

        let lf = Self::trapezoidal(lf);

        let hf = Self::trapezoidal(hf);

        let vlf = Self::trapezoidal(vlf);

        Self {
            lf: lf * T::from(2).unwrap(),
            hf: hf * T::from(2).unwrap(),
            vlf: vlf * T::from(2).unwrap(),
        }
    }

    /// Computes frequency-domain metrics (LF, HF, VLF) from raw RR intervals by first interpolating them.
    ///
    /// The `compute` function computes the frequency-domain metrics of heart rate variability (HRV) based on raw RR intervals by first interpolating them to a uniform 4 Hz rate.
    /// It utilizes the `compute_sampled` function to estimate the spectral density of the interpolated RR intervals, providing the LF, HF, and VLF values.
    ///
    /// # Parameters
    /// - `rr_intervals`: A slice of raw RR intervals, representing the time differences between successive heartbeats.
    /// - `rate`: The sampling rate used to interpolate the RR intervals, expressed in Hz.
    ///
    /// # Returns
    /// A `FrequencyMetrics` struct containing the computed LF, HF, and VLF values, which represent the power in their respective frequency bands.
    ///
    /// # Example
    /// ```
    /// use cardio_rs::frequency_domain::FrequencyMetrics;
    /// use cardio_rs::utils::test_data::RR_INTERVALS;
    ///
    /// let frequency_metrics = FrequencyMetrics::compute(RR_INTERVALS, 4.);
    /// ```
    ///
    /// # Notes
    /// This method interpolates the RR intervals using a specified rate (defaulted to 4 Hz) before passing them to the `compute_sampled` method for spectral analysis.
    pub fn compute(rr_intervals: &[T], rate: T) -> Self {
        let x = rr_intervals
            .iter()
            .scan(T::from(0).unwrap(), |s, &i| {
                *s += i;
                Some(*s)
            })
            .map(|i| (i - rr_intervals[0]) / T::from(1_000).unwrap())
            .collect::<Vec<T>>();
        let step = *x.last().unwrap() * rate;
        let xp = (0..step.to_f64().unwrap() as u64).map(|i| T::from(i).unwrap() / rate);

        let sampled_rr_intervals = Self::interpolate(&x, rr_intervals, xp);
        Self::compute_sampled(&sampled_rr_intervals, rate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_data::RR_INTERVALS;
    use approx::{AbsDiffEq, RelativeEq, UlpsEq, assert_relative_eq};

    impl<T: AbsDiffEq> AbsDiffEq for FrequencyMetrics<T>
    where
        T::Epsilon: Copy,
    {
        type Epsilon = T::Epsilon;

        fn default_epsilon() -> T::Epsilon {
            T::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
            T::abs_diff_eq(&self.lf, &other.lf, epsilon)
                && T::abs_diff_eq(&self.hf, &other.hf, epsilon)
                && T::abs_diff_eq(&self.vlf, &other.vlf, epsilon)
        }
    }

    impl<T: RelativeEq> RelativeEq for FrequencyMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_relative() -> T::Epsilon {
            T::default_max_relative()
        }

        fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
            T::relative_eq(&self.lf, &other.lf, epsilon, max_relative)
                && T::relative_eq(&self.hf, &other.hf, epsilon, max_relative)
                && T::relative_eq(&self.vlf, &other.vlf, epsilon, max_relative)
        }
    }

    impl<T: UlpsEq> UlpsEq for FrequencyMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_ulps() -> u32 {
            T::default_max_ulps()
        }

        fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
            T::ulps_eq(&self.lf, &other.lf, epsilon, max_ulps)
                && T::ulps_eq(&self.hf, &other.hf, epsilon, max_ulps)
                && T::ulps_eq(&self.vlf, &other.vlf, epsilon, max_ulps)
        }
    }

    #[test]
    fn test_frequency_metrics_with_resampling() {
        let freq_params = FrequencyMetrics::compute(RR_INTERVALS, 4.);

        assert_relative_eq!(
            FrequencyMetrics {
                lf: 3134.3763575489256,
                hf: 300.3896422597976,
                vlf: 430.1935931595116,
            },
            freq_params,
            max_relative = 0.15,
        );
    }

    #[test]
    fn test_frequency_metrics_no_resampling() {
        let freq_params = FrequencyMetrics::compute_sampled(&RR_INTERVALS.repeat(2), 4.);

        assert_relative_eq!(
            FrequencyMetrics {
                lf: 23.036923299930482,
                hf: 3676.2170268774908,
                vlf: 0.0735921788931728,
            },
            freq_params,
            max_relative = 0.03,
        );
    }
}
