//! A module providing functionality for performing window-based analysis of Heart Rate
//! Variability (HRV) metrics over non-overlapping windows of RR intervals.
//!
//! The analysis is done over fixed-duration windows, where each window is processed
//! independently, and the results for each window are stored. The user can customize the
//! pipeline for HRV metric calculation by implementing the `AnalysisPipeline` trait, or
//! they can use the default implementation that computes various HRV metrics.
//!
//! ## Key Components:
//!
//! - `WindowsAnalysisBuilder<T>`: A builder that allows you to configure and generate a
//!   `WindowsAnalysis` with user-defined data and pipeline.
//! - `AnalysisPipeline<T>`: A trait that allows the user to define custom pipelines to
//!   process RR interval data and compute various HRV metrics.
//! - `DefaultPipeline`: A default implementation of the `AnalysisPipeline` trait that computes
//!   HRV metrics using the `TimeMetrics`, `FrequencyMetrics`, and `GeometricMetrics` modules.
//! - `WindowsAnalysis<T>`: A struct representing the result of the windowed HRV analysis,
//!   containing the computed HRV metrics for each window.
//!
//! ## Example Usage:
//!
//! ### Basic Example (Using Default Pipeline):
//! ```rust
//! use cardio_rs::test_data::RR_INTERVALS;
//! use cardio_rs::windows_analysis::WindowsAnalysisBuilder;
//!
//! let rr_intervals = RR_INTERVALS.to_vec();
//! let window_size = 60_000.0; // 1-minute window size
//!
//! let windows_analysis = WindowsAnalysisBuilder::new(rr_intervals)
//!     .with_window_size(window_size)
//!     .build();
//!
//! for (i, hrv_metrics) in windows_analysis.metrics.iter().enumerate() {
//!     println!("Window {}: HRV Metrics: {:?}", i, hrv_metrics);
//! }
//! ```
//!
//! ### Custom Pipeline Example:
//! You can define your own custom analysis pipeline by implementing the `AnalysisPipeline` trait.
//! This allows you to customize how HRV metrics are computed over the RR intervals.
//!
//! ```rust
//! use cardio_rs::test_data::RR_INTERVALS;
//! use cardio_rs::{windows_analysis::{WindowsAnalysisBuilder}, HrvMetrics, geometric_domain::GeometricMetrics, processing_utils::{EctopicMethod, RRIntervals, DetectOutliers, AnalysisPipeline}, time_domain::TimeMetrics, frequency_domain::FrequencyMetrics};
//!
//! // Define a custom pipeline by implementing the AnalysisPipeline trait
//! struct CustomPipeline;
//!
//! impl AnalysisPipeline<f64> for CustomPipeline
//! {
//!     fn process(&self, data: Vec<f64>) -> HrvMetrics<f64> {
//!       let mut rr_intervals = RRIntervals::new(data);
//!       rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
//!       rr_intervals.detect_outliers(&300., &2_000.);
//!       rr_intervals.remove_outliers_ectopics();
//!
//!       let time = TimeMetrics::compute(rr_intervals.as_slice());
//!       let frequency = FrequencyMetrics::compute(rr_intervals.as_slice(), 10.);
//!       let geometric = GeometricMetrics::compute(rr_intervals.as_slice());
//!
//!       HrvMetrics {
//!           time,
//!           frequency,
//!           geometric,
//!       }
//!     }
//! }
//!
//! let rr_intervals = RR_INTERVALS.to_vec();
//!
//! // Create the window analysis with the custom pipeline and non-overlapping windows
//! let custom_pipeline = Box::new(CustomPipeline);
//! let windows_analysis = WindowsAnalysisBuilder::new(rr_intervals)
//!     .with_window_size(60_000.0) // Set the window size (1 minute)
//!     .with_pipeline(custom_pipeline) // Use the custom pipeline
//!     .build();
//!
//! // Print HRV metrics for each window
//! for (i, hrv_metrics) in windows_analysis.metrics.iter().enumerate() {
//!     println!("Window {}: HRV Metrics: {:?}", i, hrv_metrics);
//! }
//! ```

#[cfg(not(feature = "std"))]
extern crate alloc;
use crate::{HrvMetrics, processing_utils::AnalysisPipeline};
#[cfg(not(feature = "std"))]
use alloc::{boxed::Box, vec::Vec};
use core::iter::Sum;
use num::Float;

/// A struct that holds the HRV metrics for each window in the analysis.
///
/// `WindowsAnalysis<T>` stores the HRV metrics computed for each window of RR intervals.
/// The `metrics` vector contains the HRV results for all the non-overlapping windows
/// that have been processed. Each window is defined by the specified `window_size`.
pub struct WindowsAnalysis<T> {
    /// The size of the sliding window used for the analysis.
    pub window_size: T,

    /// The list of HRV metrics computed for each window.
    pub metrics: Vec<HrvMetrics<T>>,
}

/// A builder struct for configuring and constructing a `WindowsAnalysis`.
///
/// The `WindowsAnalysisBuilder` struct is used to configure the sliding window size, the analysis pipeline,
/// and the data for HRV computation. After the builder is configured, the `build()` method is used to generate
/// a `WindowsAnalysis` instance containing HRV metrics for each window.
pub struct WindowsAnalysisBuilder<T> {
    /// The RR intervals data used for the analysis.
    data: Vec<T>,

    /// The pipeline used to process the RR intervals and compute HRV metrics.
    pipeline: Box<dyn AnalysisPipeline<T>>,

    /// The size of the window used for the analysis.
    window_size: T,
}

impl<
    T: Float
        + Sum<T>
        + Copy
        + core::fmt::Debug
        + num::Signed
        + 'static
        + core::ops::AddAssign
        + core::marker::Send
        + core::marker::Sync
        + num::FromPrimitive,
> WindowsAnalysisBuilder<T>
where
    Box<dyn AnalysisPipeline<T>>: Default,
{
    /// Creates a new `WindowsAnalysisBuilder` with the provided data.
    ///
    /// # Arguments
    /// * `data` - A vector of RR intervals to be analyzed.
    ///
    /// # Returns
    /// Returns a new `WindowsAnalysisBuilder` instance.
    pub fn new(data: Vec<T>) -> Self {
        Self {
            data,
            window_size: T::from(60_000).unwrap(),
            pipeline: Default::default(),
        }
    }

    /// Sets the window size for the analysis.
    ///
    /// # Arguments
    /// * `window_size` - The size of each window, in milliseconds.
    ///
    /// # Returns
    /// Returns the builder instance for chaining.
    pub fn with_window_size(mut self, window_size: T) -> Self {
        self.window_size = window_size;
        self
    }

    /// Sets the pipeline used to process the RR intervals and compute HRV metrics.
    ///
    /// # Arguments
    /// * `pipeline` - The pipeline used to process the data.
    ///
    /// # Returns
    /// Returns the builder instance for chaining.
    pub fn with_pipeline(mut self, pipeline: Box<dyn AnalysisPipeline<T>>) -> Self {
        self.pipeline = pipeline;
        self
    }

    /// Builds the `WindowsAnalysis` struct based on the current configuration of the builder.
    ///
    /// This method processes the data into windows, computes the HRV metrics for each window,
    /// and stores the results in a `WindowsAnalysis` instance.
    ///
    /// # Returns
    /// Returns a `WindowsAnalysis` containing the HRV metrics for each window.
    pub fn build(&self) -> WindowsAnalysis<T> {
        let metrics = self
            .data
            .iter()
            .enumerate()
            .scan((T::from(0).unwrap(), 0), |s, (index, &i)| {
                s.0 += i;
                if index == self.data.len() - 1 {
                    Some(Some((s.1, self.data.len())))
                } else if s.0 <= self.window_size {
                    Some(None)
                } else {
                    s.0 = T::from(0).unwrap();
                    let prev = s.1;
                    s.1 = index;
                    Some(Some((prev, s.1)))
                }
            })
            .flatten()
            .map(|i| self.pipeline.process(self.data[i.0..i.1].to_vec()))
            .collect();

        WindowsAnalysis {
            metrics,
            window_size: self.window_size,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::RR_INTERVALS;

    #[test]
    fn test_window_size() {
        let data = RR_INTERVALS.to_vec();
        let hrv = WindowsAnalysisBuilder::new(data.clone()).build();
        assert_eq!(hrv.metrics.len(), 5);

        let hrv = WindowsAnalysisBuilder::new(data.clone())
            .with_window_size(30_000.)
            .build();
        assert_eq!(hrv.metrics.len(), 10);
    }

    #[test]
    fn test_windows_analysis() {
        let mut data = RR_INTERVALS.to_vec();
        data.extend_from_within(..);
        let hrv = WindowsAnalysisBuilder::new(data.clone())
            .with_window_size(597_000. / 2.)
            .build();
        assert_eq!(hrv.metrics[0], hrv.metrics[1]);
    }
}
