# HRV Metrics

A Rust library for computing heart rate variability (HRV) time-domain and frequency-domain metrics from RR intervals. The library provides functions for computing various HRV parameters, such as SDNN, RMSSD, LF, HF, VLF, and others, which are commonly used to assess autonomic nervous system function and heart health.

Additionally, the library offers functions for **ECG and PPG data processing and cleaning** to support accurate HRV measurements from different data sources.

## Features

- **Time-domain HRV metrics**: Compute SDNN, RMSSD, PNN50, SDSD, and mean heart rate (HR).
- **Frequency-domain HRV metrics**: Compute LF, HF, and VLF using Welch's method for power spectral density estimation.
- **ECG and PPG data processing**: Preprocess raw ECG and PPG signals for accurate HRV computation, including filtering, denoising, and R-wave peak detection for ECG and pulse detection for PPG.
- **Supports both raw and interpolated RR intervals**: Functions for both raw RR intervals and interpolated RR intervals with a specified sampling rate.
- **`no_std` compatible**: The library is designed to work in `no_std` environments, making it suitable for embedded systems.

## Time-domain Metrics

The following time-domain HRV metrics are computed:
- **RMSSD**: Root Mean Square of Successive Differences. Measures short-term HRV.
- **SDNN**: Standard Deviation of NN Intervals. Measures overall HRV.
- **PNN50**: Percentage of Successive RR Intervals > 50 ms. Indicates parasympathetic nervous system activity.
- **SDSD**: Standard Deviation of Successive Differences. Measures short-term HRV by calculating the standard deviation of differences between successive RR intervals.
- **Mean HR**: The average heart rate computed from RR intervals.
- **AVNN**: Average NN Interval, which represents the mean of RR intervals.

## Frequency-domain Metrics

The following frequency-domain HRV metrics are computed using Welch’s method for spectral density estimation:
- **LF**: Low Frequency power. Reflects the balance between sympathetic and parasympathetic nervous systems (0.04 - 0.15 Hz).
- **HF**: High Frequency power. Primarily associated with parasympathetic nervous system activity (0.15 - 0.40 Hz).
- **VLF**: Very Low Frequency power. Reflects long-term regulatory processes (0.003 - 0.04 Hz).

## ECG and PPG Data Processing and Cleaning

The library also provides functions for processing raw ECG and PPG data, including:
- **ECG Data Processing**: Includes algorithms for filtering and detecting R-wave peaks for RR interval extraction, which can be used to compute HRV metrics.
- **PPG Data Processing**: Includes preprocessing techniques for PPG signals, such as pulse detection and signal cleaning, which can also be used to extract RR intervals and calculate HRV metrics.

## Example

```rust
use cardio_rs::processing_utils::{RRIntervals, EctopicMethod, DetectOutliers};
use cardio_rs::time_domain::TimeMetrics;
use cardio_rs::frequence_domain::FrequenceMetrics;

let mut rr_intervals = RRIntervals::new(vec![800.0, 850.0, 900.0, 600.0, 800.0, 820.0, 840.0]);
rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
rr_intervals.detect_outliers(&300., &2_000.);
rr_intervals.remove_outliers_ectopics();

let frequency_metrics = FrequenceMetrics::compute(rr_intervals.as_slice());
println!("{:?}", frequency_metrics);

let time_metrics = TimeMetrics::compute(rr_intervals.as_slice());
println!("{:?}", time_metrics);
```
