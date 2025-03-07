#![cfg(feature = "std")]

use core::iter::Sum;
use interp::{InterpMode, interp_slice};
use num::Float;
use welch_sde::{Build, SpectralDensity};

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct FrequenceMetrics<T> {
    pub lf: T,
    pub hf: T,
    pub vlf: T,
}

impl<
    T: Float
        + Sum<T>
        + Copy
        + core::fmt::Debug
        + welch_sde::Signal
        + num::Signed
        + 'static
        + num::FromPrimitive,
> FrequenceMetrics<T>
{
    fn trapezoidal(x: &[T], y: &[T]) -> T {
        let mut sum = T::from(0).unwrap();
        for i in 0..(x.len() - 1) {
            let dx = x[i + 1] - x[i];
            sum += T::from(0.5).unwrap() * (y[i] + y[i + 1]) * dx;
        }
        sum
    }

    pub fn compute_sampled(sampled_rr_intervals: &[T], rate: T) -> Self {
        // TODO change welch crate to match scipy welch windows
        let n_seg = std::cmp::max(1, sampled_rr_intervals.len() / 150);
        let mean = sampled_rr_intervals.iter().copied().sum::<T>()
            / T::from_usize(sampled_rr_intervals.len()).unwrap();
        let sampled_rr_intervals: Vec<T> = sampled_rr_intervals.iter().map(|&i| i - mean).collect();
        let builder = welch_sde::Builder::new(&sampled_rr_intervals)
            .sampling_frequency(rate)
            .overlap(0.5)
            .n_segment(n_seg);

        let welch: SpectralDensity<T> = builder.build();
        let psd = welch.periodogram();

        let lf: Vec<(T, T)> = psd
            .frequency()
            .iter()
            .zip(psd.iter())
            .filter_map(|(&f, &psd_value)| {
                if f >= T::from(0.004).unwrap() && f < T::from(0.15).unwrap() {
                    Some((f, psd_value))
                } else {
                    None
                }
            })
            .collect();

        let hf: Vec<(T, T)> = psd
            .frequency()
            .iter()
            .zip(psd.iter())
            .filter_map(|(&f, &psd_value)| {
                if f >= T::from(0.15).unwrap() && f < T::from(0.4).unwrap() {
                    Some((f, psd_value))
                } else {
                    None
                }
            })
            .collect();

        let vlf: Vec<(T, T)> = psd
            .frequency()
            .iter()
            .zip(psd.iter())
            .filter_map(|(&f, &psd_value)| {
                if f >= T::from(0.003).unwrap() && f < T::from(0.04).unwrap() {
                    Some((f, psd_value))
                } else {
                    None
                }
            })
            .collect();

        let lf = Self::trapezoidal(
            &lf.iter().map(|&(f, _)| f).collect::<Vec<T>>(),
            &lf.iter()
                .map(|&(_, psd_value)| psd_value)
                .collect::<Vec<T>>(),
        );

        let hf = Self::trapezoidal(
            &hf.iter().map(|&(f, _)| f).collect::<Vec<T>>(),
            &hf.iter()
                .map(|&(_, psd_value)| psd_value)
                .collect::<Vec<T>>(),
        );

        let vlf = Self::trapezoidal(
            &vlf.iter().map(|&(f, _)| f).collect::<Vec<T>>(),
            &vlf.iter()
                .map(|&(_, psd_value)| psd_value)
                .collect::<Vec<T>>(),
        );

        Self {
            lf: lf * T::from(2).unwrap(),
            hf: hf * T::from(2).unwrap(),
            vlf: vlf * T::from(2).unwrap(),
        }
    }

    pub fn compute(rr_intervals: &[T]) -> Self {
        let x = rr_intervals
            .iter()
            .scan(T::from(0).unwrap(), |s, &i| {
                *s += i;
                Some(*s)
            })
            .map(|i| (i - rr_intervals[0]) / T::from(1_000).unwrap())
            .collect::<Vec<T>>();
        let rate = T::from(4).unwrap();
        let step = *x.last().unwrap() * rate;
        let xp = (0..step.to_f64().unwrap() as u64)
            .map(|i| T::from(i).unwrap() / rate)
            .collect::<Vec<T>>();

        let sampled_rr_intervals = interp_slice(&x, &rr_intervals, &xp, &InterpMode::Extrapolate);
        Self::compute_sampled(&sampled_rr_intervals, rate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_data::test_data::RR_INTERVALS;
    use approx::AbsDiffEq;
    use approx::RelativeEq;
    use approx::UlpsEq;
    use approx::assert_relative_eq;

    impl<T: AbsDiffEq> AbsDiffEq for FrequenceMetrics<T>
    where
        T::Epsilon: Copy,
    {
        type Epsilon = T::Epsilon;

        fn default_epsilon() -> T::Epsilon {
            T::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
            T::abs_diff_eq(&self.lf, &other.lf, epsilon)
                && T::abs_diff_eq(&self.hf, &other.hf, epsilon)
                && T::abs_diff_eq(&self.vlf, &other.vlf, epsilon)
        }
    }

    impl<T: RelativeEq> RelativeEq for FrequenceMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_relative() -> T::Epsilon {
            T::default_max_relative()
        }

        fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
            T::relative_eq(&self.lf, &other.lf, epsilon, max_relative)
                && T::relative_eq(&self.hf, &other.hf, epsilon, max_relative)
                && T::relative_eq(&self.vlf, &other.vlf, epsilon, max_relative)
        }
    }

    impl<T: UlpsEq> UlpsEq for FrequenceMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_ulps() -> u32 {
            T::default_max_ulps()
        }

        fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
            T::ulps_eq(&self.lf, &other.lf, epsilon, max_ulps)
                && T::ulps_eq(&self.hf, &other.hf, epsilon, max_ulps)
                && T::ulps_eq(&self.vlf, &other.vlf, epsilon, max_ulps)
        }
    }

    #[test]
    fn test_frequency_metrics() {
        let freq_params = FrequenceMetrics::compute(&RR_INTERVALS);

        // TODO: comparaison is done https://github.com/Aura-healthcare/hrv-analysis
        // seems to not be the same as neurokit2
        assert_relative_eq!(
            FrequenceMetrics {
                lf: 3134.3763575489256,
                hf: 300.3896422597976,
                vlf: 430.1935931595116,
            },
            freq_params,
            epsilon = 500.,
        );
    }
}
