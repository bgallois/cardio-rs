# ❤️⚡ **CARDIO-RS** ⚡❤️  
*A Rust library for Heart Rate Variability (HRV) analysis*  

[![Docs.rs](https://docs.rs/cardio-rs/badge.svg)](https://docs.rs/cardio-rs)  [![CI](https://github.com/bgallois/cardio-rs/actions/workflows/test.yml/badge.svg)](https://github.com/bgallois/cardio-rs/actions/) [![Crates.io](https://img.shields.io/crates/v/cardio-rs.svg)](https://crates.io/crates/cardio-rs)

📊 **Compute HRV metrics** in time-domain, frequency-domain, and non-linear domain from RR intervals with ease!  
💓 **Preprocess ECG & PPG raw data** for accurate HRV analysis!  
📟 **Compatible with embedded systems** using `no_std`!

---

## ✨ **Key Features**  

- ✅ **Time-domain HRV metrics**: SDNN, RMSSD, PNN50, SDSD, CVSD, and more!  
- ✅ **Frequency-domain HRV metrics**: LF, HF, VLF (using Welch’s method).  
- ✅ **Geometric HRV metrics**: Triangular Index & TINN.  
- ✅ **Non-linear HRV metrics**: SampEn, DFA, LZC, and more!  
- ✅ **ECG & PPG preprocessing**: Filtering, denoising, peak detection.  
- ✅ **`no_std` compatibility** for embedded systems.  
- ✅ **Window-based analysis**: Split data into windows for segment-based HRV calculations.  
- ✅ **Live HRV analysis**: Real-time processing of RR intervals with customizable pipelines.

---

## 📏 **Time-domain HRV Metrics**  

- 🔹 **RMSSD** – Short-term HRV measurement.  
- 🔹 **SDNN** – Overall HRV measurement.  
- 🔹 **PNN50** – Percentage of successive RR intervals > 50ms.  
- 🔹 **SDSD** – Standard deviation of successive RR intervals.  
- 🔹 **Mean HR** – Average heart rate.  
- 🔹 **AVNN** – Mean of all RR intervals.  
- 🔹 **CVSD** – RMSSD divided by the mean RR interval.

---

## 📐 **Geometric HRV Metrics**  

- 📊 **Triangular Index** – Measures RR interval distribution.  
- 📊 **TINN** (Triangular Interpolation of NN Interval Histogram) – Estimates RR variability.

---

## 🎵 **Frequency-domain HRV Metrics**  

**Computed using Welch’s method** for spectral density estimation:

- 🎼 **LF (Low Frequency)** – Sympathetic and parasympathetic balance (**0.04 - 0.15 Hz**).  
- 🎼 **HF (High Frequency)** – Parasympathetic activity (**0.15 - 0.40 Hz**).  
- 🎼 **VLF (Very Low Frequency)** – Long-term regulatory processes (**0.003 - 0.04 Hz**).  

---

## 🧩 **Non-linear HRV Metrics**  

Non-linear HRV metrics provide insights into the complexity and self-organization of the heart's autonomic regulation. They are useful for capturing dynamics that cannot be fully understood through linear analysis alone.

- 🔹 **SampEn (Sample Entropy)** – Measures the regularity and unpredictability of a time series.  
- 🔹 **DFA (Detrended Fluctuation Analysis)** – Evaluates long-range correlations and scaling behavior in the RR intervals.  
- 🔹 **LZC (Lempel-Ziv Complexity)** – Quantifies the complexity of a binary sequence derived from RR intervals.

---

## 🏥 **ECG & PPG Data Processing**  

- 📡 **ECG Processing**: R-wave peak detection, filtering, and RR extraction.  
- 💡 **PPG Processing**: Pulse detection & signal cleaning for HRV computation.

---

## 🪟 **Window Analysis**  

Cardio-rs supports **window analysis**, which allows you to split data into segments of a defined size. This is useful for performing HRV analysis on specific time periods or segments of data.

### Example:

```rust
use cardio_rs::{ windows_analysis::WindowsAnalysisBuilder };

let data = vec![800., 810., 790., 765., 780., 800., 810., 795., 770., 785.];

// Create the window analysis builder
let analysis = WindowsAnalysisBuilder::new(data)
    .with_window_size(3_000.)
    .build();

// Get HRV metrics for each window
println!("{:?}", analysis.metrics);
```

---

## 🔄 **Live Analysis**  

The **live analysis** module allows for real-time HRV calculations by dynamically processing incoming RR intervals. You can define a custom analysis pipeline or use the default pipeline for standard HRV metrics.

### Example:

```rust
use cardio_rs::{ live_analysis::TimeQueue, HrvMetrics };
use cardio_rs::utils::test_data::RR_INTERVALS;

let mut queue = TimeQueue::new(60_000);  // Set the time window (in milliseconds)
let rr_intervals = RR_INTERVALS.to_vec();

for interval in rr_intervals {
    queue.push(interval);
}

let hrv = queue.get_hrv();
println!("Calculated HRV: {:?}", hrv);
```

---

## 🦀 **Example Usage**  

### 🚀 **Compute HRV Metrics & Load ECG Data**  

```rust
use cardio_rs::{ standard_analysis, HrvMetrics };

let path = "tests/ecg.csv";
let signal = "ECG_Raw";
let rate = 1000.;

let hrv_metrics = standard_analysis!(path, signal, rate);
println!("{:?}", hrv_metrics);
```

```rust
#[cfg(feature = "std")] {
use cardio_rs::{
    processing_utils::{ RRIntervals, EctopicMethod, DetectOutliers },
    time_domain::TimeMetrics,
    non_linear::NonLinearMetrics,
    geometric_domain::GeometricMetrics,
    frequency_domain::FrequencyMetrics,
    io_utils::{ DataBuilder },
    standard_analysis,
    HrvMetrics,
};

let path = "tests/ecg.csv";
let signal = "ECG_Raw";
let time = "Time";
let data = DataBuilder::new(path.into(), signal.into())
    .with_time(time.into())
    .build()
    .unwrap();

let hrv_metrics = standard_analysis!(data.clone().get_rr());
println!("{:?}", hrv_metrics);

// Manual analysis for more control
let mut rr_intervals = RRIntervals::new(data.get_rr());

rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
rr_intervals.detect_outliers(&300., &2_000.);
rr_intervals.remove_outliers_ectopics();

let time_metrics = TimeMetrics::compute(rr_intervals.as_slice());
println!("{:?}", time_metrics);

let frequency_metrics = FrequencyMetrics::compute(rr_intervals.as_slice(), 10.);
println!("{:?}", frequency_metrics);

let geo_metrics = GeometricMetrics::compute(rr_intervals.as_slice());
println!("{:?}", geo_metrics);

let non_linear_metrics = NonLinearMetrics::compute_default(&rr_intervals);
println!("{:?}", non_linear_metrics);
}
```

---

🔥 **Start analyzing HRV today with cardio-rs!** 🚀💓
