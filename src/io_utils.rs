//! A module providing functionality for processing and analyzing time series data from CSV files, particularly for extracting RR intervals
//! from physiological signals such as ECG and PPG.
//!
//! It defines the `Data` struct for handling the processed data and the
//! `DataBuilder` struct for building and configuring data sources.
//!
//! ## Key Components:
//! - `Data`: The main struct representing the processed time series data. It provides methods for collecting the data
//!   and calculating RR intervals based on detected peaks.
//! - `DataBuilder`: A builder struct that helps construct a `Data` object, allowing users to specify signal columns,
//!   time columns, or a sampling rate for the data processing.
//!
//! ## Features:
//! - **RR Interval Extraction**: The `Data` struct provides a method `get_rr` that calculates RR intervals based on detected
//!   peaks in the signal, using a time series of `time` and `signal` columns or a custom rate for time generation.
//! - **Flexible Data Construction**: The `DataBuilder` struct provides methods to configure the dataset, allowing users to
//!   specify time columns, rate, or signal columns, and ensuring the data is correctly loaded and processed.
//! - **Peak Detection**: The module uses the `find_peaks` library for detecting peaks in the signal, which are essential for
//!   calculating the RR intervals.
//!
//! ## Example Usage:
//!
//! ```rust
//! use cardio_rs::io_utils::{DataBuilder, Data};
//!
//! let path = "tests/ecg.csv";
//! let signal = "ECG_Raw";
//! let time = "Time";
//!
//! let data = DataBuilder::new(path.into(), signal.into())
//!     .with_time(time.into())
//!     .build()
//!     .unwrap();
//!
//! let rr_intervals = data.get_rr();
//! ```
//!
//! ```rust
//! use cardio_rs::io_utils::{DataBuilder, Data};
//!
//! let path = "tests/ecg.csv";
//! let signal = "ECG_Raw";
//!
//! let data = DataBuilder::new(path.into(), signal.into())
//!     .with_rate(1_000.)
//!     .build()
//!     .unwrap();
//!
//! let rr_intervals = data.get_rr();
//! ```
#![cfg(feature = "std")]
use find_peaks::PeakFinder;
use polars::{
    datatypes::DataType,
    frame::DataFrame,
    prelude::{CsvReadOptions, NamedFrom, PolarsError, SerReader},
    series::Series,
};
use polars_lazy::{
    frame::IntoLazy,
    prelude::{LazyFrame, col, lit},
};
use std::path::PathBuf;

/// `Data` represents a collection of time series data processed lazily.
///
/// This struct wraps a `LazyFrame`, which is a deferred computation
/// representation of a dataset. It provides functionality to perform
/// operations like collecting the data into a `DataFrame` or extracting
/// RR intervals from physiological signals.
#[derive(Clone)]
pub struct Data {
    /// Raw data.
    data: LazyFrame,
}

impl Data {
    /// Collects the lazy data into an eager `DataFrame`.
    ///
    /// This method forces the execution of the lazy computation, returning
    /// the resulting `DataFrame`.
    ///
    /// # Returns
    /// - A `DataFrame` containing the collected data.
    pub fn collect(self) -> DataFrame {
        self.data.collect().unwrap()
    }

    /// Extracts RR intervals from the time series data.
    ///
    /// This method uses peak detection
    /// to compute the RR intervals. The RR intervals are returned as a vector
    /// of floating-point numbers in milliseconds.
    ///
    /// # Returns
    /// - A `Vec<f64>` representing the RR intervals in milliseconds.
    pub fn get_rr(self) -> Vec<f64> {
        let df = self.data.clone().collect().unwrap();

        let x: Vec<f64> = df
            .column("time")
            .unwrap()
            .cast(&DataType::Float64)
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect();

        let signal = df
            .column("signal")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        let ps = PeakFinder::new_with_x(&signal, &x)
            .with_min_prominence(0.75)
            .with_min_distance(0.33)
            .find_peaks();

        let mut ps: Vec<f64> = ps.iter().map(|i| x[i.middle_position()]).collect();
        ps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        ps.windows(2).map(|i| (i[1] - i[0]) * 1000.0).collect()
    }
}

/// Builder pattern for constructing `Data` objects.
///
/// The `DataBuilder` allows for the specification of input data files,
/// signal column names, and options for time or rate. It ensures that
/// the data is correctly loaded, pre-processed, and formatted for analysis.
pub struct DataBuilder {
    file: PathBuf,
    signal: String,
    time: Option<String>,
    rate: Option<f64>,
}

impl DataBuilder {
    /// Creates a new `DataBuilder` instance.
    ///
    /// # Arguments
    /// - `file`: The path to the CSV file containing the data.
    /// - `signal`: The label of the column containing the signal data.
    ///
    /// # Returns
    /// - A new `DataBuilder` instance.
    pub fn new(file: PathBuf, signal: String) -> Self {
        DataBuilder {
            file,
            signal,
            time: None,
            rate: None,
        }
    }

    /// Sets the time column for the `DataBuilder`.
    ///
    /// This method allows for specifying a column that represents time.
    ///
    /// # Arguments
    /// - `time`: The label of the time column in the dataset.
    ///
    /// # Panics
    /// Panics if the `rate` is already set.
    ///
    /// # Returns
    /// - The updated `DataBuilder` instance.
    pub fn with_time(mut self, time: String) -> Self {
        if self.rate.is_some() {
            panic!("Rate is already set");
        }
        self.time = Some(time);
        self
    }

    /// Sets the rate for the `DataBuilder`.
    ///
    /// This method allows specifying a sampling rate to generate a time column
    /// when a time column is not provided.
    ///
    /// # Arguments
    /// - `rate`: The sampling rate in Hz.
    ///
    /// # Panics
    /// Panics if the `time` is already set.
    ///
    /// # Returns
    /// - The updated `DataBuilder` instance.
    pub fn with_rate(mut self, rate: f64) -> Self {
        if self.time.is_some() {
            panic!("Time is already set");
        }
        self.rate = Some(rate);
        self
    }

    /// Builds a `Data` object based on the current configuration.
    ///
    /// This method reads the CSV file, selects the signal column, optionally
    /// adds the time column or generates one based on the rate, and returns
    /// the constructed `Data` object.
    ///
    /// # Returns
    /// - A `Result<Data, PolarsError>` representing the result of the build process.
    pub fn build(self) -> Result<Data, PolarsError> {
        let mut data = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(self.file))?
            .finish()?
            .lazy();
        let mut selected_cols = vec![
            ((col(&self.signal) - col(&self.signal).mean()) / col(&self.signal).std(1))
                .alias("signal"),
        ];

        if let Some(time) = self.time {
            selected_cols.push(col(time).alias("time"));
        } else if let Some(rate) = self.rate {
            let len_expr = data
                .clone()
                .select([col(&self.signal)])
                .collect()
                .unwrap()
                .height();
            let time_series: Vec<f64> = (0..len_expr).map(|i| i as f64 / rate).collect();
            let time_series = Series::new("time".into(), time_series);
            let time_series_lit = lit(time_series).alias("time");
            selected_cols.push(time_series_lit);
        }

        data = data.select(&selected_cols);
        Ok(Data { data })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_process_ppg() {
        let path = "tests/ppg.csv";
        let signal = "PPG_Raw";
        let time = "Time";
        let data = DataBuilder::new(path.into(), signal.into())
            .with_time(time.into())
            .build()
            .unwrap();
        let rr = data.get_rr();

        let reference = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(path.into()))
            .unwrap()
            .finish()
            .unwrap()
            .lazy();
        let peaks = reference
            .filter(col("PPG_Peaks").eq(lit(1)))
            .select([col("Time").cast(DataType::Float64)])
            .collect()
            .unwrap()
            .column("Time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        let reference_rr: Vec<f64> = peaks.windows(2).map(|w| (w[1] - w[0]) * 1_000.).collect();

        assert_eq!(rr, reference_rr);
    }

    #[test]
    fn test_process_ecg() {
        let path = "tests/ecg.csv";
        let signal = "ECG_Raw";
        let time = "Time";
        let data = DataBuilder::new(path.into(), signal.into())
            .with_time(time.into())
            .build()
            .unwrap();
        let rr = data.get_rr();

        let reference = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(path.into()))
            .unwrap()
            .finish()
            .unwrap()
            .lazy();
        let peaks = reference
            .filter(col("ECG_Peaks").eq(lit(1)))
            .select([col("Time").cast(DataType::Float64)])
            .collect()
            .unwrap()
            .column("Time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        let reference_rr: Vec<f64> = peaks.windows(2).map(|w| (w[1] - w[0]) * 1_000.).collect();

        assert_eq!(rr, reference_rr);
    }

    #[test]
    fn test_process_rate() {
        let path = "tests/ecg.csv";
        let signal = "ECG_Raw";
        let data = DataBuilder::new(path.into(), signal.into())
            .with_rate(1_000.)
            .build()
            .unwrap();
        let rr = data.get_rr();

        let reference = CsvReadOptions::default()
            .with_has_header(true)
            .try_into_reader_with_file_path(Some(path.into()))
            .unwrap()
            .finish()
            .unwrap()
            .lazy();
        let peaks = reference
            .filter(col("ECG_Peaks").eq(lit(1)))
            .select([col("Time").cast(DataType::Float64)])
            .collect()
            .unwrap()
            .column("Time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        let reference_rr: Vec<f64> = peaks.windows(2).map(|w| (w[1] - w[0]) * 1_000.).collect();

        assert_eq!(rr, reference_rr);
    }
}
