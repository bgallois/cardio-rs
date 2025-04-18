//! A module providing functionality for detecting and handling outliers and ectopic beats in a series of RR intervals.
//!
//! It defines the `RRIntervals` struct which contains a vector of RR intervals, and optional fields to store indices of
//! detected outliers and ectopic beats. The module also provides a trait `DetectOutliers` which can be implemented to detect
//! outliers and ectopic beats using custom methods, while default implementations for detecting outliers and ectopic beats
//! using common algorithms are also provided. The module includes methods for removing detected outliers and ectopics from
//! the RR intervals data.
//!
//! ## Key Components:
//! - `RRIntervals`: The main struct representing a collection of RR intervals and optional outlier/ectopic detections.
//! - `DetectOutliers`: A trait for detecting outliers and ectopics with customizable implementations.
//! - `EctopicMethod`: An enum that defines different methods for detecting ectopic beats in RR intervals.
//!
//! ## Features:
//! - **Detecting Outliers**: Outliers in the RR intervals are detected by comparing each interval to a provided range (`lowest_rr`, `highest_rr`).
//! - **Detecting Ectopics**: Ectopic beats are detected using methods like the Karlsson method, which compares the mean of adjacent intervals.
//! - **Flexible Customization**: Users can implement their own methods for detecting outliers and ectopics by implementing the `DetectOutliers` trait.
//! - **Removing Outliers and Ectopics**: Detected outliers and ectopics can be removed from the data using the `remove_outliers_ectopics` method, leaving only the valid RR intervals.
//!
//! ## Example Usage:
//!
//! ```rust
//! use cardio_rs::processing_utils::{RRIntervals, EctopicMethod, DetectOutliers};
//!
//! let mut rr_intervals = RRIntervals::new(vec![800.0, 850.0, 3000.0, 600.0, 800.0]);
//! rr_intervals.detect_outliers(&300.0, &2000.0);  // Detect outliers based on specified range
//! rr_intervals.detect_ectopics(EctopicMethod::Karlsson);  // Detect ectopic beats using the Karlsson method
//! rr_intervals.remove_outliers_ectopics();  // Remove outliers and ectopics from the RR intervals data
//! ```
//!
//! ## Trait and Struct Documentation:
//! - `RRIntervals<T>`: A struct representing a sequence of RR intervals along with optional detected outliers and ectopics.
//! - `DetectOutliers<T>`: A trait for detecting outliers and ectopics. Custom implementations can be provided by the user.
//! - `EctopicMethod`: An enum for specifying different methods to detect ectopic beats (currently only `Karlsson` is supported).
//!

#[cfg(not(feature = "std"))]
extern crate alloc;
#[cfg(not(feature = "std"))]
use alloc::{boxed::Box, vec, vec::Vec};
use core::{
    iter::Sum,
    ops::{Deref, DerefMut},
};
use num::Float;

/// Enum representing different methods for detecting ectopic beats in RR intervals.
///
/// This enum provides various algorithms for identifying and removing ectopic beats
/// based on predefined statistical criteria.
pub enum EctopicMethod {
    /// Karlsson method for detecting ectopic beats in RR intervals.
    ///
    /// An RR interval is considered ectopic if it differs by more than 20%
    /// from the mean of the previous and next RR intervals.
    Karlsson,

    /// Acar method for detecting ectopic beats in RR intervals.
    ///
    /// An RR interval is considered ectopic if it differs by more than 20%
    /// from the mean of the last 9 RR intervals.
    Acar,
}

/// Struct representing RR intervals and associated outlier and ectopic detection results.
///
/// This struct contains a vector of RR intervals and optional fields to store indices of
/// detected outliers and ectopic beats. It provides methods to detect these events
/// and remove them from the data as needed.
#[derive(Debug)]
pub struct RRIntervals<T> {
    /// A vector of RR intervals representing the data.
    rr_intervals: Vec<T>,

    /// An optional vector storing indices of detected outliers in the RR intervals.
    outliers: Option<Vec<usize>>,

    /// An optional vector storing indices of detected ectopic beats in the RR intervals.
    ectopics: Option<Vec<usize>>,
}

/// Trait for detecting outliers and ectopics in RR intervals.
///
/// This trait allows custom implementations of methods for detecting outliers and ectopic beats in RR intervals.
/// Any type that implements this trait can provide different detection methods based on specific requirements.
pub trait DetectOutliers<T> {
    fn detect_outliers(&mut self, lowest_rr: &T, highest_rr: &T);
    fn detect_ectopics(&mut self, method: EctopicMethod);
}

impl<T: Float + Sum<T> + Copy + core::fmt::Debug + num::Signed + 'static + num::FromPrimitive>
    DetectOutliers<T> for RRIntervals<T>
{
    /// Detects outliers in the RR intervals based on the given `lowest_rr` and `highest_rr`.
    ///
    /// The method identifies RR intervals that are either lower than `lowest_rr` or higher than `highest_rr`
    /// as outliers. These outliers are stored in the `outliers` field of the struct.
    ///
    /// # Arguments
    ///
    /// * `lowest_rr` - The minimum acceptable RR interval.
    /// * `highest_rr` - The maximum acceptable RR interval.
    ///
    /// # Example
    ///
    /// ```rust
    /// use cardio_rs::processing_utils::{RRIntervals, DetectOutliers};
    ///
    /// let mut rr_intervals = RRIntervals::new(vec![800.0, 850.0, 3000.0, 600.0, 800.0]);
    /// rr_intervals.detect_outliers(&300.0, &2000.0);
    /// ```
    fn detect_outliers(&mut self, lowest_rr: &T, highest_rr: &T) {
        let outliers: Vec<usize> = self
            .iter()
            .enumerate()
            .filter_map(|(index, value)| {
                if lowest_rr > value || value > highest_rr {
                    Some(index)
                } else {
                    None
                }
            })
            .collect();
        self.outliers = if outliers.is_empty() {
            None
        } else {
            Some(outliers)
        };
    }

    /// Detects ectopic beats in the RR intervals using a specified method.
    ///
    /// This method detects ectopic beats based on the given `method` (e.g., the `Karlsson` method),
    /// and stores the detected indices of ectopics in the `ectopics` field of the struct.
    ///
    /// # Arguments
    ///
    /// * `method` - The method used to detect ectopic beats.
    ///
    /// # Example
    ///
    /// ```rust
    /// use cardio_rs::processing_utils::{RRIntervals, EctopicMethod, DetectOutliers};
    ///
    /// let mut rr_intervals = RRIntervals::new(vec![800.0, 850.0, 900.0, 600.0, 800.0]);
    /// rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
    /// ```
    fn detect_ectopics(&mut self, method: EctopicMethod) {
        let ectopics: Vec<usize> = match method {
            EctopicMethod::Karlsson => (0..self.len() - 2)
                .filter_map(|i| {
                    let mean = (self[i] + self[i + 2]) / T::from(2).unwrap();
                    if (mean - self[i + 1]).abs() >= T::from(0.2).unwrap() * mean {
                        Some(i + 1)
                    } else {
                        None
                    }
                })
                .collect(),
            EctopicMethod::Acar => (9..self.len())
                .filter(|&i| {
                    let mean = (self[i - 9..i].iter().cloned().sum::<T>()) / T::from(9).unwrap();
                    (mean - self[i]).abs() >= T::from(0.2).unwrap() * mean
                })
                .collect(),
        };
        self.ectopics = if ectopics.is_empty() {
            None
        } else {
            Some(ectopics)
        };
    }
}

impl<T> Deref for RRIntervals<T> {
    type Target = Vec<T>;

    /// Deref implementation for accessing the underlying `Vec<T>` of RR intervals.
    ///
    /// This implementation allows easy read-only access to the vector of RR intervals.
    fn deref(&self) -> &Self::Target {
        &self.rr_intervals
    }
}

impl<T> DerefMut for RRIntervals<T> {
    /// DerefMut implementation for mutable access to the underlying `Vec<T>` of RR intervals.
    ///
    /// This implementation allows modifying the vector of RR intervals directly.
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.rr_intervals
    }
}

impl<T: Float + Sum<T> + Copy + core::fmt::Debug + num::Signed + 'static + num::FromPrimitive>
    RRIntervals<T>
{
    /// Creates a new instance of `RRIntervals` from a vector of RR intervals.
    ///
    /// # Arguments
    ///
    /// * `rr_intervals` - A vector of RR intervals.
    ///
    /// # Returns
    ///
    /// A new `RRIntervals` instance containing the provided RR intervals, with `None` values for outliers and ectopics.
    pub fn new(rr_intervals: Vec<T>) -> Self {
        Self {
            rr_intervals,
            outliers: None,
            ectopics: None,
        }
    }

    /// Removes outliers and ectopics from the RR intervals.
    ///
    /// This method removes any elements from the RR intervals vector at indices that
    /// correspond to detected outliers or ectopics. After removal, the `outliers` and `ectopics`
    /// fields are reset to `None`.
    pub fn remove_outliers_ectopics(&mut self) {
        self.rr_intervals = self
            .iter()
            .enumerate()
            .filter_map(|(i, j)| {
                if self.outliers.as_ref().unwrap_or(&vec![]).contains(&i)
                    || self.ectopics.as_ref().unwrap_or(&vec![]).contains(&i)
                {
                    None
                } else {
                    Some(*j)
                }
            })
            .collect::<Vec<T>>();
        self.ectopics = None;
        self.outliers = None;
    }
}

/// A trait for processing HRV data through an analysis pipeline.
pub trait AnalysisPipeline<T> {
    /// Processes a vector of RR intervals to compute HRV metrics.
    ///
    /// # Arguments
    /// * `data` - A vector of RR intervals representing the time between heartbeats for a given window.
    ///
    /// # Returns
    /// Returns a struct containing the computed HRV metrics (time-domain, frequency-domain, and geometric).
    fn process(&self, data: Vec<T>) -> crate::HrvMetrics<T>;
}

/// A default implementation of the `AnalysisPipeline` trait that computes HRV metrics using a predefined method.
struct DefaultPipeline();
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
        + Into<f64>
        + num::FromPrimitive,
> AnalysisPipeline<T> for DefaultPipeline
{
    fn process(&self, data: Vec<T>) -> crate::HrvMetrics<T> {
        let mut rr_intervals = RRIntervals::new(data);
        rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
        rr_intervals.detect_outliers(&T::from(300).unwrap(), &T::from(2_000).unwrap());
        rr_intervals.remove_outliers_ectopics();

        let time = crate::time_domain::TimeMetrics::compute(rr_intervals.as_slice());
        let frequency = crate::frequency_domain::FrequencyMetrics::compute(
            rr_intervals.as_slice(),
            T::from(4).unwrap(),
        );
        let geometric = crate::geometric_domain::GeometricMetrics::compute(rr_intervals.as_slice());
        let non_linear =
            crate::non_linear::NonLinearMetrics::compute_default(rr_intervals.as_slice());

        crate::HrvMetrics {
            time,
            frequency,
            geometric,
            non_linear,
        }
    }
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
> Default for Box<dyn AnalysisPipeline<T>>
where
    f64: From<T>,
{
    /// Returns a default instance of the `DefaultPipeline`.
    ///
    /// This method returns the default `DefaultPipeline` that computes HRV metrics using default methods.
    fn default() -> Self {
        Box::new(DefaultPipeline())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::test_data::RR_INTERVALS;

    #[test]
    fn test_remove_none() {
        let mut rr_intervals = RRIntervals::new(RR_INTERVALS.to_vec());
        rr_intervals.remove_outliers_ectopics();
        assert_eq!(RR_INTERVALS, *rr_intervals);
    }

    #[test]
    fn test_remove_ectopics_karlsson() {
        let mut rr_intervals = RRIntervals::new(vec![800., 850., 900., 600., 800., 820., 840.]);
        rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
        assert_eq!(Some(vec![2, 3]), rr_intervals.ectopics);
        rr_intervals.remove_outliers_ectopics();
        assert_eq!(vec![800., 850., 800., 820., 840.], *rr_intervals);
    }

    #[test]
    fn test_remove_ectopics_acar() {
        let mut rr_intervals = RRIntervals::new(vec![
            800., 850., 3000., 600., 800., 820., 240., 800., 850., 3000., 600., 800., 820., 240.,
        ]);
        rr_intervals.detect_ectopics(EctopicMethod::Acar);
        assert_eq!(Some(vec![9, 10, 11, 13]), rr_intervals.ectopics);
        rr_intervals.remove_outliers_ectopics();
        assert_eq!(
            vec![800., 850., 3000., 600., 800., 820., 240., 800., 850., 820.],
            *rr_intervals
        );
    }

    #[test]
    fn test_remove_outliers() {
        let mut rr_intervals = RRIntervals::new(vec![800., 850., 3000., 600., 800., 820., 240.]);
        rr_intervals.detect_outliers(&300., &2_000.);
        assert_eq!(Some(vec![2, 6]), rr_intervals.outliers);
        rr_intervals.remove_outliers_ectopics();
        assert_eq!(vec![800., 850., 600., 800., 820.], *rr_intervals);
    }

    #[test]
    fn test_remove_combined() {
        let mut rr_intervals = RRIntervals::new(vec![800., 850., 3000., 600., 800., 820., 240.]);
        rr_intervals.outliers = Some(vec![2, 3, 4]);
        rr_intervals.ectopics = Some(vec![4, 5]);
        rr_intervals.remove_outliers_ectopics();
        assert_eq!(vec![800., 850., 240.], *rr_intervals);
    }
}
