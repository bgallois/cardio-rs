use approx::assert_relative_eq;
use cardio_rs::{
    frequency_domain::FrequencyMetrics,
    geometric_domain::GeometricMetrics,
    io_utils::DataBuilder,
    processing_utils::{DetectOutliers, EctopicMethod, RRIntervals},
    time_domain::TimeMetrics,
};

#[test]
fn test_real_ecg() {
    let path = "tests/mimic_perform_af_csv/mimic_perform_af_001_data.csv";
    let signal = "ECG";
    let time = "Time";
    let data = DataBuilder::new(path.into(), signal.into())
        .with_time(time.into())
        .build()
        .unwrap();
    let mut rr_intervals = RRIntervals::new(data.get_rr());

    rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
    rr_intervals.detect_outliers(&300., &2_000.);
    rr_intervals.remove_outliers_ectopics();

    let time_metrics = TimeMetrics::compute(rr_intervals.as_slice());

    let frequency_metrics = FrequencyMetrics::compute(rr_intervals.as_slice(), 4.);

    let geo_metrics = GeometricMetrics::compute(rr_intervals.as_slice());

    let precision = 0.05;
    assert_relative_eq!(
        time_metrics.rmssd,
        78.4745076238388,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.sdnn,
        70.2880299250242,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.pnn50,
        0.49375459221160,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.mean_hr,
        93.44226608425353,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.sdsd,
        78.47450212026048,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.cvsd,
        0.12070874146616935,
        max_relative = precision
    );

    assert_relative_eq!(
        geo_metrics.triangular_index,
        14.48936170212766,
        max_relative = precision
    );

    let precision = 0.2;
    assert_relative_eq!(
        frequency_metrics.hf,
        1459.8731219921774,
        max_relative = precision
    );
    assert_relative_eq!(
        frequency_metrics.lf,
        1278.4757430583918,
        max_relative = precision
    );
    assert_relative_eq!(
        frequency_metrics.vlf,
        300.76349399523303,
        max_relative = precision
    );
}

#[test]
fn test_real_ppg() {
    let path = "tests/mimic_perform_af_csv/mimic_perform_af_001_data.csv";
    let signal = "PPG";
    let time = "Time";
    let data = DataBuilder::new(path.into(), signal.into())
        .with_time(time.into())
        .build()
        .unwrap();
    let mut rr_intervals = RRIntervals::new(data.get_rr());

    rr_intervals.detect_ectopics(EctopicMethod::Karlsson);
    rr_intervals.detect_outliers(&300., &2_000.);
    rr_intervals.remove_outliers_ectopics();

    let time_metrics = TimeMetrics::compute(rr_intervals.as_slice());

    let frequency_metrics = FrequencyMetrics::compute(rr_intervals.as_slice(), 4.);

    let geo_metrics = GeometricMetrics::compute(rr_intervals.as_slice());

    let precision = 0.05;
    assert_relative_eq!(
        time_metrics.rmssd,
        81.35094813338031,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.sdnn,
        70.93915105482384,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.pnn50,
        0.4903175832687839,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.mean_hr,
        92.12599294303625,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.sdsd,
        81.35094600926843,
        max_relative = precision
    );
    assert_relative_eq!(
        time_metrics.cvsd,
        0.12346635332381918,
        max_relative = precision
    );

    let precision = 0.2;

    assert_relative_eq!(geo_metrics.triangular_index, 15.2, max_relative = precision);
    assert_relative_eq!(
        frequency_metrics.hf,
        1618.31099240089,
        max_relative = precision
    );
    assert_relative_eq!(
        frequency_metrics.lf,
        1377.9704183067502,
        max_relative = precision
    );
    assert_relative_eq!(
        frequency_metrics.vlf,
        387.16583492213647,
        max_relative = precision
    );
}
