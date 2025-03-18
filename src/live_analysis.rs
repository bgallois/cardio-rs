//! A module providing the `TimeQueue` struct, which manages a fixed size time series of RR intervals
//! and allows for real-time processing of heart rate variability (HRV) metrics. It utilizes
//! a sliding window of RR intervals and processes them using a user-defined or default
//! pipeline for HRV calculation.
//!
//! # Examples
//!
//! Basic usage with the default pipeline:
//! ```
//! use cardio_rs::test_data::RR_INTERVALS;
//! use cardio_rs::live_analysis::TimeQueue;
//! let mut queue = TimeQueue::new(60_000);
//! let rr_intervals = RR_INTERVALS.to_vec();
//! for interval in rr_intervals {
//!     queue.push(interval);
//! }
//! let hrv = queue.get_hrv();
//! println!("Calculated HRV: {:?}", hrv);
//! ```
//!
//! Custom pipeline usage:
//! ```
//! use cardio_rs::{live_analysis::TimeQueue, windows_analysis::AnalysisPipeline, HrvMetrics};
//! struct CustomPipeline;
//! impl AnalysisPipeline<f64> for CustomPipeline {
//!     fn process(&self, data: Vec<f64>) -> HrvMetrics<f64> {
//!         HrvMetrics::default()
//!     }
//! }
//! use cardio_rs::test_data::RR_INTERVALS;
//! let mut queue = TimeQueue::new(60_000);
//! queue.set_pipeline(Box::new(CustomPipeline));
//! let rr_intervals = RR_INTERVALS.to_vec();
//! for interval in rr_intervals {
//!     queue.push(interval);
//! }
//! let hrv = queue.get_hrv();
//! println!("Calculated HRV with custom pipeline: {:?}", hrv);
//! ```
use crate::windows_analysis::AnalysisPipeline;
use num::Float;
#[cfg(feature = "std")]
use std::collections::VecDeque;
#[cfg(not(feature = "std"))]
extern crate alloc;
#[cfg(not(feature = "std"))]
use alloc::{boxed::Box, collections::vec_deque::VecDeque};

/// A struct to manage a sliding window of RR intervals and process heart rate variability (HRV) metrics.
/// It stores a collection of RR intervals in `data` and calculates HRV using a custom or default pipeline.
pub struct TimeQueue<T> {
    /// A `VecDeque` that stores the current RR intervals in the sliding window.
    data: VecDeque<T>,

    /// The total window time (in ms). Once this time is exceeded, the oldest RR interval is removed.
    time: T,

    /// The accumulated time that tracks how much time (in ms) has passed in the current window.
    current_time: T,

    /// The pipeline used to process the RR intervals and compute HRV metrics. It can be set to a custom pipeline.
    pipeline: Box<dyn AnalysisPipeline<T>>,
}

impl<
    T: Float
        + core::iter::Sum<T>
        + Copy
        + 'static
        + core::fmt::Debug
        + num::Signed
        + core::ops::AddAssign
        + core::marker::Send
        + core::marker::Sync
        + core::ops::SubAssign
        + Into<f64>
        + num::FromPrimitive,
> TimeQueue<T>
where
    Box<dyn AnalysisPipeline<T>>: Default,
{
    /// Creates a new `TimeQueue` with a sliding window of size `time`. A default pipeline is used.
    ///
    /// # Arguments
    /// * `time` - The length of the time window (in ms).
    pub fn new(time: usize) -> Self {
        Self {
            data: VecDeque::with_capacity(4 * time), // At most 250 bpm
            time: T::from(time).unwrap(),
            current_time: T::from(0).unwrap(),
            pipeline: Default::default(),
        }
    }

    /// Sets a custom pipeline to be used for HRV calculation.
    ///
    /// # Arguments
    /// * `pipeline` - A boxed `AnalysisPipeline` that implements the HRV processing logic.
    pub fn set_pipeline(&mut self, pipeline: Box<dyn AnalysisPipeline<T>>) {
        self.pipeline = pipeline;
    }

    /// Pushes a new RR interval to the queue. The queue maintains the sliding window of RR intervals.
    /// Once the current time exceeds the defined window time, the oldest interval is removed.
    pub fn push(&mut self, rr_interval: T) {
        self.current_time += rr_interval;
        self.data.push_back(rr_interval);
        if self.current_time >= self.time {
            if let Some(deleted) = self.data.pop_front() {
                self.current_time -= deleted;
            }
        }
    }

    /// Returns the current slice of the stored RR intervals.
    pub fn get(&self) -> &[T] {
        self.data.as_slices().0
    }

    /// Processes the current RR intervals using the pipeline and returns the HRV metrics.
    ///
    /// This function calls the `process` method on the custom or default pipeline to calculate HRV.
    pub fn get_hrv(&self) -> super::HrvMetrics<T> {
        self.pipeline.process(self.get().to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        frequency_domain::FrequencyMetrics, geometric_domain::GeometricMetrics,
        processing_utils::RRIntervals, test_data::RR_INTERVALS, time_domain::TimeMetrics,
    };
    #[cfg(not(feature = "std"))]
    use alloc::vec::Vec;

    #[test]
    fn test_queue() {
        let mut data = RR_INTERVALS.to_vec();
        data.extend_from_within(..);
        let mut queue = TimeQueue::new(298_500usize);
        let mut win_0: crate::HrvMetrics<f64> = Default::default();
        for (i, j) in data.iter().enumerate() {
            queue.push(*j);
            if i == data.len() / 2 - 1 {
                win_0 = queue.get_hrv();
            }
        }

        let win_1 = queue.get_hrv();
        assert_eq!(win_0, win_1);
    }
    #[test]
    fn test_queue_custom_pipeline() {
        let mut data = RR_INTERVALS.to_vec();
        data.extend_from_within(..);
        let mut queue = TimeQueue::new(298_500usize);

        struct Pipeline();
        impl AnalysisPipeline<f64> for Pipeline {
            fn process(&self, data: Vec<f64>) -> crate::HrvMetrics<f64> {
                let rr_intervals = RRIntervals::new(data);

                let time = TimeMetrics::compute(rr_intervals.as_slice());
                let frequency = FrequencyMetrics::compute(rr_intervals.as_slice(), 4.);
                let geometric = GeometricMetrics::compute(rr_intervals.as_slice());

                crate::HrvMetrics {
                    time,
                    frequency,
                    geometric,
                }
            }
        }
        queue.set_pipeline(Box::new(Pipeline()));

        let mut win_0: crate::HrvMetrics<f64> = Default::default();
        for (i, j) in data.iter().enumerate() {
            queue.push(*j);
            if i == data.len() / 2 - 1 {
                win_0 = queue.get_hrv();
            }
        }

        let win_1 = queue.get_hrv();
        assert_eq!(win_0, win_1);
    }
}
