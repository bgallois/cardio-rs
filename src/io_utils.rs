#![cfg(feature = "std")]
use find_peaks::PeakFinder;
use polars::frame::DataFrame;
use polars::prelude::CsvReadOptions;
use polars::prelude::PolarsError;
use polars::prelude::SerReader;
use polars::series::ops::NullBehavior;
use polars_lazy::frame::IntoLazy;
use polars_lazy::prelude::LazyFrame;
use polars_lazy::prelude::col;
use std::path::PathBuf;

pub struct Data {
    data: LazyFrame,
}

impl Data {
    pub fn collect(self) -> DataFrame {
        self.data.collect().unwrap()
    }

    pub fn get_rr(self) -> Vec<f64> {
        let df = self.data.clone().collect().unwrap();

        let rate = self
            .data
            .clone()
            .select([col("time").diff(1, NullBehavior::Ignore).mean()])
            .collect()
            .unwrap()
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .get(0)
            .unwrap();

        let x = df
            .column("time")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        let signal = df
            .column("signal")
            .unwrap()
            .f64()
            .unwrap()
            .into_no_null_iter()
            .collect::<Vec<f64>>();

        let ps = PeakFinder::new_with_x(&signal, &x)
            .with_min_height(0.2)
            .with_min_distance(rate)
            .find_peaks();

        let mut ps: Vec<f64> = ps.iter().map(|i| x[i.middle_position()]).collect();
        ps.sort_by(|a, b| a.partial_cmp(b).unwrap());

        ps.windows(2).map(|i| (i[1] - i[0]) * 1000.0).collect()
    }
}

pub struct DataBuilder {
    file: PathBuf,
    signal: String,
    time: Option<String>,
}

impl DataBuilder {
    pub fn new(file: PathBuf, signal: String) -> Self {
        DataBuilder {
            file,
            signal,
            time: None,
        }
    }

    pub fn with_time(mut self, time: String) -> Self {
        self.time = Some(time);
        self
    }

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
        }
        data = data.select(&selected_cols);
        Ok(Data { data })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read() {
        let path = "tests/mimic_perform_af_csv/mimic_perform_af_001_data.csv";
        let signal = "PPG";
        let time = "Time";
        let data = DataBuilder::new(path.into(), signal.into())
            .with_time(time.into())
            .build()
            .unwrap();
        let _ = data.get_rr();
    }
}
