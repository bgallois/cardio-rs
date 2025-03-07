extern crate alloc;

use alloc::vec::Vec;
use core::iter::Sum;
use num::Float;

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct TimeMetrics<T> {
    pub rmssd: T,
    pub sdnn: T,
    pub pnn50: T,
}

impl<T: Float + Sum<T> + Copy + core::fmt::Debug> TimeMetrics<T> {
    pub fn compute(rr_intervals: &[T]) -> Self {
        let rr_diffs: Vec<T> = rr_intervals.windows(2).map(|i| i[1] - i[0]).collect();
        let rr_intervals_mean: T =
            rr_intervals.iter().copied().sum::<T>() / T::from(rr_intervals.len()).unwrap();

        let variance = rr_intervals
            .iter()
            .map(|&x| (x - rr_intervals_mean) * (x - rr_intervals_mean))
            .sum::<T>()
            / T::from(rr_diffs.len()).unwrap();

        let sdnn = variance.sqrt();

        let rmssd =
            (rr_diffs.iter().map(|&x| x * x).sum::<T>() / T::from(rr_diffs.len()).unwrap()).sqrt();

        let pnn50 = T::from(
            rr_diffs
                .iter()
                .filter(|&&x| x.abs() > T::from(50.0).unwrap())
                .count(),
        )
        .unwrap()
            / T::from(rr_diffs.len()).unwrap();
        Self { sdnn, rmssd, pnn50 }
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
    impl<T: AbsDiffEq> AbsDiffEq for TimeMetrics<T>
    where
        T::Epsilon: Copy,
    {
        type Epsilon = T::Epsilon;

        fn default_epsilon() -> T::Epsilon {
            T::default_epsilon()
        }

        fn abs_diff_eq(&self, other: &Self, epsilon: T::Epsilon) -> bool {
            T::abs_diff_eq(&self.rmssd, &other.rmssd, epsilon)
                && T::abs_diff_eq(&self.sdnn, &other.sdnn, epsilon)
                && T::abs_diff_eq(&self.pnn50, &other.pnn50, epsilon)
        }
    }

    impl<T: RelativeEq> RelativeEq for TimeMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_relative() -> T::Epsilon {
            T::default_max_relative()
        }

        fn relative_eq(&self, other: &Self, epsilon: T::Epsilon, max_relative: T::Epsilon) -> bool {
            T::relative_eq(&self.rmssd, &other.rmssd, epsilon, max_relative)
                && T::relative_eq(&self.sdnn, &other.sdnn, epsilon, max_relative)
                && T::relative_eq(&self.pnn50, &other.pnn50, epsilon, max_relative)
        }
    }

    impl<T: UlpsEq> UlpsEq for TimeMetrics<T>
    where
        T::Epsilon: Copy,
    {
        fn default_max_ulps() -> u32 {
            T::default_max_ulps()
        }

        fn ulps_eq(&self, other: &Self, epsilon: T::Epsilon, max_ulps: u32) -> bool {
            T::ulps_eq(&self.rmssd, &other.rmssd, epsilon, max_ulps)
                && T::ulps_eq(&self.sdnn, &other.sdnn, epsilon, max_ulps)
                && T::ulps_eq(&self.pnn50, &other.pnn50, epsilon, max_ulps)
        }
    }

    #[test]
    fn test_metrics() {
        let time_params = TimeMetrics::compute(&RR_INTERVALS);

        assert_relative_eq!(
            TimeMetrics {
                rmssd: 59.99807873965272,
                sdnn: 69.54843274772418,
                pnn50: 0.414985590778098,
            },
            time_params,
            epsilon = 1e-14
        );
    }
}
