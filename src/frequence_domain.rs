#![cfg(feature = "std")]

use core::iter::Sum;
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
    pub fn compute(sampled_rr_intervals: &[T], rate: T) -> Self {
        let welch: SpectralDensity<T> =
            SpectralDensity::builder(&sampled_rr_intervals, rate).build();
        let psd = welch.periodogram();
        let delta_f = T::from(rate).unwrap() / T::from(psd.frequency().len()).unwrap();

        let vlf = psd
            .iter()
            .zip(psd.frequency())
            .filter_map(|(&p, f)| {
                if f >= T::from(0.0033).unwrap() && f < T::from(0.04).unwrap() {
                    return Some(p * delta_f);
                } else {
                    None
                }
            })
            .sum::<T>();

        let lf = psd
            .iter()
            .zip(psd.frequency())
            .filter_map(|(&p, f)| {
                if f >= T::from(0.04).unwrap() && f < T::from(0.15).unwrap() {
                    return Some(p * delta_f);
                } else {
                    None
                }
            })
            .sum::<T>();

        let hf = psd
            .iter()
            .zip(psd.frequency())
            .filter_map(|(&p, f)| {
                if f >= T::from(0.15).unwrap() && f < T::from(0.4).unwrap() {
                    return Some(p * delta_f);
                } else {
                    None
                }
            })
            .sum::<T>();

        Self { lf, hf, vlf }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::AbsDiffEq;
    use approx::RelativeEq;
    use approx::UlpsEq;
    use approx::assert_relative_eq;
    use rand_distr::{Distribution, Normal};

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
        let rr_dist = Normal::new(800., 50.).unwrap();
        let rr_intervals: Vec<_> = (0..1_000_000)
            .map(|_| rr_dist.sample(&mut rand::rng()))
            .collect();
        let rate = 1.0;
        let freq_params = FrequenceMetrics::compute(&rr_intervals, rate);

        assert_relative_eq!(
            FrequenceMetrics {
                lf: 549.,
                hf: 1252.,
                vlf: 184.,
            },
            freq_params,
            epsilon = 10.
        );
    }

    #[test]
    fn test_frequency_metrics_rate() {
        let rr_dist = Normal::new(800., 50.).unwrap();
        let rr_intervals: Vec<_> = (0..10_000_000)
            .map(|_| rr_dist.sample(&mut rand::rng()))
            .collect();
        let rate = 10.0;
        let freq_params = FrequenceMetrics::compute(&rr_intervals, rate);

        assert_relative_eq!(
            FrequenceMetrics {
                lf: 54.9,
                hf: 125.2,
                vlf: 18.4,
            },
            freq_params,
            epsilon = 1.
        );
    }
}
