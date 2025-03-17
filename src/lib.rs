#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]

pub mod frequency_domain;
pub mod geometric_domain;
pub mod io_utils;
pub mod processing_utils;
pub mod test_data;
pub mod time_domain;
pub mod welch;

/// A struct that contains the Heart Rate Variability (HRV) metrics,
/// divided into time-domain, frequency-domain, and geometric-domain metrics.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct HrvMetrics<T> {
    /// Time-domain HRV metrics (`TimeMetrics<T>`).
    pub time: crate::time_domain::TimeMetrics<T>,
    /// Frequency-domain HRV metrics (`FrequencyMetrics<T>`).
    pub frequency: crate::frequency_domain::FrequencyMetrics<T>,
    /// Geometric-domain HRV metrics (`GeometricMetrics<T>`).
    pub geometric: crate::geometric_domain::GeometricMetrics<T>,
}

/// A macro that performs a standard HRV analysis on a given input, returning the HRV metrics
/// in time-domain, frequency-domain, and geometric-domain categories.
///
/// This macro computes HRV metrics using a series of transformations:
/// 1. It processes RR intervals by detecting ectopic beats and outliers.
/// 2. Computes time-domain metrics, frequency-domain metrics, and geometric-domain metrics.
///
/// The macro accepts two forms:
///
/// - **Form 1**: Takes a reference to an object containing RR intervals (e.g., `RRIntervals`).
///
/// ```rust
/// use cardio_rs::standard_analysis;
/// use cardio_rs::HrvMetrics;
///
/// let rr_intervals = vec![800., 810., 780., 850., 900.];
/// let hrv_metrics = standard_analysis!(rr_intervals);
/// ```
///
/// - **Form 2**: Takes a file path, signal, and sample rate to build a `Data` object and then runs the analysis on that data.
///
/// ```rust
/// #[cfg(feature = "std")]
/// {
/// use cardio_rs::standard_analysis;
/// use cardio_rs::HrvMetrics;
///
/// let path = "tests/ecg.csv";
/// let signal = "ECG_Raw";
/// let rate = 1000.;
/// let hrv_metrics = cardio_rs::standard_analysis!(path, signal, rate);
/// }
/// ```
/// Where `path` is the file path to the data, `signal` is the name of the signal, and `rate`
/// is the sample rate of the signal.
///
/// The macro relies on the following modules automatically imported for computing HRV metrics:
/// - `cardio_rs::frequency_domain::FrequencyMetrics`
/// - `cardio_rs::geometric_domain::GeometricMetrics`
/// - `cardio_rs::time_domain::TimeMetrics`
/// - `cardio_rs::processing_utils::{DetectOutliers, EctopicMethod, RRIntervals}`
/// - `cardio_rs::io_utils::{Data, DataBuilder}`
#[macro_export]
macro_rules! standard_analysis {
    ( $x:expr ) => {{
        use cardio_rs::{
            frequency_domain::FrequencyMetrics,
            geometric_domain::GeometricMetrics,
            processing_utils::{DetectOutliers, EctopicMethod, RRIntervals},
            time_domain::TimeMetrics,
        };
        let mut rr_intervals = RRIntervals::new($x);
        rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
        rr_intervals.detect_outliers(&300., &2_000.);
        rr_intervals.remove_outliers_ectopics();

        let time = TimeMetrics::compute(rr_intervals.as_slice());
        let frequency = FrequencyMetrics::compute(rr_intervals.as_slice(), 10.);
        let geometric = GeometricMetrics::compute(rr_intervals.as_slice());

        HrvMetrics {
            time,
            frequency,
            geometric,
        }
    }};
    ( $path:expr, $signal:expr, $rate:expr ) => {{
        #[cfg(feature = "std")]
        {
            use cardio_rs::io_utils::{Data, DataBuilder};
            let data = DataBuilder::new($path.into(), $signal.into())
                .with_rate($rate)
                .build()
                .unwrap();
            standard_analysis!(data.get_rr())
        }
    }};
}
