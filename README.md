# ❤️⚡ CARDIO-RS ⚡❤️  
*A Rust library for Heart Rate Variability (HRV) analysis!*  

[![Docs.rs](https://docs.rs/cardio-rs/badge.svg)](https://docs.rs/cardio-rs)  [![CI](https://github.com/bgallois/cardio-rs/actions/workflows/test.yml/badge.svg)](https://github.com/bgallois/cardio-rs/actions/)[![Crates.io](https://img.shields.io/crates/v/cardio-rs.svg)](https://crates.io/crates/cardio-rs)

📊 **Compute HRV time-domain & frequency-domain & geometric-domain metrics** from RR intervals with ease!  
💓 **Preprocess ECG & PPG raw data** for accurate HRV analysis!  
📟 **Supports embedded systems** with `no_std` (work in progress 🚧)!  

> **Note:** Not all functionalities are `no_std` yet, but I’m working on it! 🛠️  

---

## ✨ Features  

✅ **Time-domain HRV metrics**: SDNN, RMSSD, PNN50, SDSD, CVSD, and more!  
✅ **Frequency-domain HRV metrics**: LF, HF, VLF using Welch’s method!  
✅ **Geometric-domain HRV metrics**: Triangular Index & TINN!  
✅ **ECG & PPG preprocessing**: Filtering, denoising, peak detection!  
✅ **Raw & interpolated RR intervals** supported!  
✅ **🚀 `no_std` compatibility** (in progress)!  

---

## 📏 Time-domain HRV Metrics  

- 🔹 **RMSSD** – Measures short-term HRV.  
- 🔹 **SDNN** – Measures overall HRV.  
- 🔹 **PNN50** – Percentage of successive RR intervals > 50ms.  
- 🔹 **SDSD** – Standard deviation of successive differences.  
- 🔹 **Mean HR** – Average heart rate.  
- 🔹 **AVNN** – Mean of all RR intervals.  
- 🔹 **CVSD** – RMSSD divided by the mean RR interval.  

---

## 📐 Geometric-domain HRV Metrics  

- 📊 **Triangular Index** – Measures RR interval distribution.  
- 📊 **TINN (Triangular Interpolation of NN Interval Histogram)** – Estimates RR variability.  

---

## 🎵 Frequency-domain HRV Metrics  

📡 Computed using **Welch’s method** for spectral density estimation:  

- 🎼 **LF (Low Frequency)** – Sympathetic & parasympathetic balance (**0.04 - 0.15 Hz**).  
- 🎼 **HF (High Frequency)** – Parasympathetic activity (**0.15 - 0.40 Hz**).  
- 🎼 **VLF (Very Low Frequency)** – Long-term regulatory processes (**0.003 - 0.04 Hz**).  

---

## 🏥 ECG & PPG Data Processing  

📡 **ECG Processing** – R-wave peak detection, filtering, and RR extraction!  
💡 **PPG Processing** – Pulse detection & signal cleaning for HRV computation!  

---

## 🦀 Example Usage  

### 🚀 Compute HRV Metrics & Load ECG Data  

```rust
use cardio_rs::{
    processing_utils::{RRIntervals, EctopicMethod, DetectOutliers},
    time_domain::TimeMetrics,
    geometric_domain::GeometricMetrics,
    frequency_domain::FrequencyMetrics,
    io_utils::{DataBuilder, Data},
};

let path = "tests/ecg.csv";
let signal = "ECG_Raw";
let time = "Time";
let data = DataBuilder::new(path.into(), signal.into())
    .with_time(time.into())
    .build()
    .unwrap();
let mut rr_intervals = RRIntervals::new(data.get_rr());

rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
rr_intervals.detect_outliers(&300., &2_000.);
rr_intervals.remove_outliers_ectopics();

let freq_metrics = FrequencyMetrics::compute(rr_intervals.as_slice());
println!("{:?}", freq_metrics);

let time_metrics = TimeMetrics::compute(rr_intervals.as_slice());
println!("{:?}", time_metrics);

let geo_metrics = GeometricMetrics::compute(rr_intervals.as_slice());
println!("{:?}", geo_metrics);
```

🔥 **Start analyzing HRV today with cardio-rs!** 🚀💓
