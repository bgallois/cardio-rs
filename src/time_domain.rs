use arrayvec::ArrayVec;
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
        let rr_diffs: ArrayVec<T, 100> = rr_intervals.windows(2).map(|i| i[1] - i[0]).collect();
        let rr_diffs_mean: T =
            rr_diffs.iter().copied().sum::<T>() / T::from(rr_diffs.len()).unwrap();

        let variance = rr_diffs
            .iter()
            .map(|&x| (x - rr_diffs_mean) * (x - rr_diffs_mean))
            .sum::<T>()
            / T::from(rr_intervals.len() - 1usize).unwrap();

        let sdnn = variance.sqrt();
        let rmssd = (rr_diffs.iter().map(|&x| x * x).sum::<T>()
            / T::from(rr_intervals.len() - 1usize).unwrap())
        .sqrt();
        let pnn50 = T::from(
            rr_diffs
                .iter()
                .filter(|&&x| x.abs() > T::from(50.0).unwrap())
                .count(),
        )
        .unwrap()
            / T::from(rr_intervals.len() - 1usize).unwrap();
        Self { sdnn, rmssd, pnn50 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
        let rr_intervals: &[f64] = &[800., 810., 790., 850., 795., 900., 860., 950.];
        let time_params = TimeMetrics::compute(rr_intervals);

        assert_relative_eq!(
            TimeMetrics {
                rmssd: 63.075918f64,
                sdnn: 59.324428f64,
                pnn50: 0.57142857f64,
            },
            time_params,
            epsilon = 1e-7
        );
    }
}
