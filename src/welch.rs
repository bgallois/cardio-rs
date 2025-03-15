//! This module provides an implementation of spectral analysis tools,
//! including the Hann window and Welch's method for power spectral density estimation.
//!
//! Unlike existing Rust implementations, this approach closely follows SciPy’s methodology,
//! ensuring compatibility and familiar behavior for users migrating from Python.
//!
//! This implementation should (and will) eventually be extracted into its own crate,
//! with additional configuration options for greater flexibility and usability.
//!
//! The module includes:
//! - `Hann`: A Hann window function for reducing spectral leakage.
//! - `Welch`: Welch's method for estimating the power spectral density.
//! - `HannBuilder` and `WelchBuilder` for constructing instances with configurable parameters.
//!
//! Future improvements should focus on extensibility, performance optimizations,
//! and additional options to match SciPy's functionality.

#[cfg(not(feature = "std"))]
extern crate alloc;
#[cfg(not(feature = "std"))]
use alloc::vec;
#[cfg(not(feature = "std"))]
use alloc::vec::Vec;
use core::{cmp, f64::consts::PI};
use num::{Complex, Float};
#[cfg(feature = "std")]
use rustfft::{Fft, FftDirection};

/// A trait representing a windowing function that can be applied to a signal.
///
/// Window functions are commonly used in signal processing to mitigate spectral leakage
/// when performing Fourier transforms. A window function modifies a signal by applying
/// a predefined weighting curve before analysis.
pub trait Window<T> {
    /// Applies the window function to the given signal.
    ///
    /// # Parameters
    /// * `signal` - A vector representing the input signal.
    ///
    /// # Returns
    /// A new iterator where each sample has been multiplied by the corresponding window weight.
    ///
    /// # Panics
    /// Implementations may panic if `signal.len()` does not match the expected window size.
    fn apply(&self, chunk: &[T]) -> impl Iterator<Item = T>;
    /// Computes the power of the window function.
    ///
    /// The power is typically the sum of squared window coefficients, which can be used
    /// for normalization purposes in signal processing.
    ///
    /// # Returns
    /// The computed power of the window function.
    fn power(&self) -> T;
    /// Computes the sum of the window function.
    ///
    /// # Returns
    /// The computed sum of the window function.
    fn sum(&self) -> T;
}

impl<T: Float + Copy + core::fmt::Debug + core::iter::Sum> Window<T> for Hann<T> {
    fn apply(&self, chunk: &[T]) -> impl Iterator<Item = T> {
        if chunk.len() != self.weights.len() {
            panic!("Signal and Window should have the same size");
        }
        let chunk_mean = chunk.iter().copied().sum::<T>() / T::from(chunk.len()).unwrap();
        chunk
            .iter()
            .zip(self.weights.iter())
            .map(move |(&a, &b)| (a - chunk_mean) * b)
    }

    fn power(&self) -> T {
        self.weights.iter().map(|&w| w * w).sum()
    }

    fn sum(&self) -> T {
        self.weights.iter().copied().sum()
    }
}

/// Represents a Hann window function.
///
/// The Hann window is commonly used in spectral analysis and signal processing
/// to reduce spectral leakage by applying a tapering function to the signal.
pub struct Hann<T> {
    weights: Vec<T>,
}

/// Builder for creating a Hann window.
///
/// This struct allows for flexible construction of a Hann window with a given size.
pub struct HannBuilder {
    n: usize,
}

impl HannBuilder {
    /// Creates a new HannBuilder with the given window size.
    ///
    /// # Arguments
    /// * `n` - The number of points in the window.
    ///
    /// # Returns
    /// A new `HannBuilder` instance.
    pub fn new(n: usize) -> Self {
        HannBuilder { n }
    }

    /// Constructs a Hann window of the specified type.
    ///
    /// # Type Parameters
    /// * `T` - A floating-point type that supports FFT operations.
    ///
    /// # Returns
    /// A `Hann<T>` instance containing the computed window weights.
    pub fn build<T>(&self) -> Hann<T>
    where
        T: Float + Copy + core::fmt::Debug + core::iter::Sum,
    {
        let weights = (0..self.n)
            .map(|i| {
                T::from(0.5 * (1.0 - ((2.0 * PI * i as f64) / (self.n as f64 - 1.0)).cos()))
                    .unwrap()
            })
            .collect::<Vec<T>>();
        Hann { weights }
    }
}

/// Represents different normalization methods for spectral analysis.
///
/// This enum is used to specify how the periodogram should be normalized when computing
/// the power spectral density (PSD) or the power spectrum.
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Normalization<T> {
    /// Normalization by sampling frequency, producing a power spectral density (PSD).
    ///
    /// When used, the resulting spectrum has units of `V²/Hz` if the input signal is in V
    /// and the sampling frequency is in Hz. The total power of the signal is obtained
    /// by integrating over all frequencies.
    Density,

    /// Normalization by the sum of squared window coefficients, producing a power spectrum.
    ///
    /// When used, the resulting spectrum has units of `V²` if the input signal is in volts.
    /// This represents the total power distributed across frequency bins.
    Spectrum,

    /// Custom normalization factor.
    ///
    /// Allows the user to specify a custom normalization value, useful for advanced applications
    /// where neither `Density` nor `Spectrum` provide the desired scaling.
    Custom(T),
}

/// A trait for computing the periodogram and corresponding frequencies of a signal.
///
/// The periodogram estimates the power spectral density (PSD) of a signal.
/// Implementations of this trait should provide methods to compute both the periodogram
/// and its associated frequency values.
pub trait Periodogram<T> {
    /// Returns an iterator over the periodogram values.
    ///
    /// The periodogram represents the power spectral density (PSD) of a signal,
    /// computed using a spectral estimation method (e.g., Welch's method).
    ///
    /// # Returns
    /// An iterator that yields the power spectral density values as `T`.
    fn periodogram(&self) -> impl Iterator<Item = T> + '_;

    /// Returns an iterator over the frequency values corresponding to the periodogram.
    ///
    /// The frequencies are typically derived based on the sampling rate and windowing
    /// used in the spectral estimation.
    ///
    /// # Returns
    /// An iterator that yields the frequency values as `T`.
    fn frequencies(&self) -> impl Iterator<Item = T> + '_;
}

impl<T: Float + Copy + core::fmt::Debug + core::iter::Sum> Periodogram<T> for Welch<T> {
    fn periodogram(&self) -> impl Iterator<Item = T> + '_ {
        self.periodogram.iter().copied()
    }

    fn frequencies(&self) -> impl Iterator<Item = T> + '_ {
        self.frequencies.iter().copied()
    }
}

/// Represents the result of Welch's method for power spectral density estimation.
///
/// Welch's method reduces noise in power spectral estimates by averaging
/// multiple periodograms computed from overlapping segments.
///
/// # Type Parameters
/// * `T` - A numeric type that supports floating-point operations.
pub struct Welch<T> {
    periodogram: Vec<T>,
    frequencies: Vec<T>,
}

/// Builder for constructing a Welch power spectral density estimator.
///
/// This struct provides a flexible way to configure the parameters for Welch's method,
/// such as segment size, overlap, and FFT size.
///
/// # Type Parameters
/// * `T` - A floating-point type that supports FFT operations.
pub struct WelchBuilder<T> {
    normalization: Normalization<T>,
    segment_size: usize,
    dft_size: usize,
    overlap_size: usize,
    fs: T,
    signal: Vec<T>,
}

impl<
    T: Float
        + Copy
        + core::fmt::Debug
        + core::marker::Sync
        + core::marker::Send
        + core::iter::Sum
        + core::ops::AddAssign
        + num::Signed
        + num::FromPrimitive
        + 'static,
> WelchBuilder<T>
{
    /// Creates a new WelchBuilder with a given input signal.
    ///
    /// # Arguments
    /// * `signal` - The input signal as a vector of type `T`.
    ///
    /// # Returns
    /// A new `WelchBuilder` instance with default parameters.
    pub fn new(signal: Vec<T>) -> Self {
        WelchBuilder {
            normalization: Normalization::Density,
            signal,
            segment_size: 256,
            dft_size: 4096,
            overlap_size: 128,
            fs: T::from(4).unwrap(),
        }
    }

    /// Sets the segment size for Welch's method.
    pub fn with_segment_size(mut self, n: usize) -> Self {
        self.segment_size = n;
        self
    }

    /// Sets the FFT size.
    #[cfg(feature = "std")]
    pub fn with_dft_size(mut self, n: usize) -> Self {
        self.dft_size = n;
        self
    }

    /// Sets the overlap size between segments.
    pub fn with_overlap_size(mut self, n: usize) -> Self {
        self.overlap_size = n;
        self
    }

    /// Sets the sampling frequency.
    pub fn with_fs(mut self, n: T) -> Self {
        self.fs = n;
        self
    }

    /// Sets the normalization.
    pub fn with_normalization(mut self, norma: Normalization<T>) -> Self {
        self.normalization = norma;
        self
    }

    /// Constructs the Welch power spectral density estimator.
    ///
    /// # Returns
    /// A `Welch<T>` instance containing the computed periodogram and frequencies.
    pub fn build(&self) -> Welch<T> {
        let window = HannBuilder::new(self.segment_size).build::<T>();
        let mut periodogram: Vec<T> = vec![T::from(0).unwrap(); self.dft_size / 2];
        let mut i = 0;
        let mut n_segments = 0;
        while i + self.segment_size < self.signal.len() {
            let chunk = &self.signal[i..cmp::min(i + self.segment_size, self.signal.len())];
            let chunk = window.apply(chunk);

            let mut buffer = vec![
                Complex {
                    re: T::from(0).unwrap(),
                    im: T::from(0).unwrap(),
                };
                self.dft_size
            ];
            chunk.enumerate().for_each(|(i, j)| {
                buffer[i].re = j;
            });

            #[cfg(feature = "std")]
            {
                let fft = rustfft::algorithm::Radix4::new(self.dft_size, FftDirection::Forward);
                fft.process(&mut buffer);
            }

            #[cfg(not(feature = "std"))]
            {
                naive_fft(&mut buffer);
            }

            let pdg: Vec<T> = buffer
                .into_iter()
                .take(self.dft_size / 2)
                .map(|i| i.norm_sqr())
                .collect();
            periodogram
                .iter_mut()
                .zip(pdg.into_iter())
                .for_each(|(a, b)| *a += b);

            i += self.segment_size - self.overlap_size;
            n_segments += 1;
        }
        let frequencies: Vec<T> = (0..self.dft_size / 2)
            .map(|i| T::from(i).unwrap() * self.fs / T::from(self.dft_size).unwrap())
            .collect();
        let norma = match self.normalization {
            Normalization::Density => {
                (window.power() * T::from(n_segments).unwrap() * self.fs).recip()
            }
            Normalization::Spectrum => {
                (window.sum() * window.sum() * T::from(n_segments).unwrap()).recip()
            }
            Normalization::Custom(e) => e * T::from(n_segments).unwrap(),
        };
        let periodogram = periodogram.into_iter().map(|p| p * norma).collect();

        Welch {
            periodogram,
            frequencies,
        }
    }
}

#[cfg(not(feature = "std"))]
fn naive_fft<T: Float + Copy + core::fmt::Debug + core::iter::Sum>(input: &mut [Complex<T>]) {
    let n = input.len();
    if n <= 1 {
        return;
    }

    let mut even: Vec<Complex<T>> = input.iter().copied().step_by(2).collect();
    let mut odd: Vec<Complex<T>> = input.iter().copied().skip(1).step_by(2).collect();

    naive_fft(&mut even);
    naive_fft(&mut odd);

    for k in 0..n / 2 {
        let twiddle = Complex::from_polar(
            T::one(),
            -T::from(2.0 * PI).unwrap() * T::from(k).unwrap() / T::from(n).unwrap(),
        ) * odd[k];
        input[k] = even[k] + twiddle;
        input[k + n / 2] = even[k] - twiddle;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::rng;
    use rand_distr::{Distribution, Normal};

    #[test]
    fn test_normal_density() {
        let normal = Normal::new(0., 1.).unwrap();
        let mut rng = rng();
        let samples: Vec<f64> = (0..100_000).map(|_| normal.sample(&mut rng)).collect();

        let welch = WelchBuilder::new(samples)
            .with_fs(4.)
            .with_overlap_size(128)
            .with_segment_size(256)
            .with_normalization(Normalization::Density)
            .build();
        assert_relative_eq!(
            welch.periodogram().sum::<f64>()
                / welch.periodogram().collect::<Vec<f64>>().len() as f64,
            0.249,
            epsilon = 1e-2
        );
    }

    #[test]
    fn test_normal_spectrum() {
        let normal = Normal::new(0., 1.).unwrap();
        let mut rng = rng();
        let samples: Vec<f64> = (0..100_000).map(|_| normal.sample(&mut rng)).collect();

        let welch = WelchBuilder::new(samples)
            .with_fs(4.)
            .with_overlap_size(128)
            .with_segment_size(256)
            .with_normalization(Normalization::Spectrum)
            .build();
        assert_relative_eq!(
            welch.periodogram().sum::<f64>()
                / welch.periodogram().collect::<Vec<f64>>().len() as f64,
            0.00585,
            epsilon = 1e-4
        );
    }

    #[test]
    fn test_frequence() {
        let normal = Normal::new(0., 1.).unwrap();
        let mut rng = rng();
        let samples: Vec<f64> = (0..100_000).map(|_| normal.sample(&mut rng)).collect();

        let welch = WelchBuilder::new(samples)
            .with_fs(400.)
            .with_overlap_size(128)
            .with_segment_size(256)
            .with_normalization(Normalization::Density)
            .build();
        assert_relative_eq!(
            welch.periodogram().sum::<f64>()
                / welch.periodogram().collect::<Vec<f64>>().len() as f64,
            0.0024846134053086084,
            epsilon = 1e-4
        );
    }

    #[cfg(feature = "std")]
    #[test]
    fn test_dtf() {
        let normal = Normal::new(0., 1.).unwrap();
        let mut rng = rng();
        let samples: Vec<f64> = (0..100_000).map(|_| normal.sample(&mut rng)).collect();

        let welch = WelchBuilder::new(samples)
            .with_fs(400.)
            .with_dft_size(512)
            .with_overlap_size(128)
            .with_segment_size(256)
            .with_normalization(Normalization::Density)
            .build();
        assert_relative_eq!(
            welch.periodogram().sum::<f64>()
                / welch.periodogram().collect::<Vec<f64>>().len() as f64,
            0.0025083335651038905,
            epsilon = 1e-4
        );
    }

    #[cfg(feature = "std")]
    #[test]
    fn test_short() {
        let normal = Normal::new(0., 1.).unwrap();
        let mut rng = rng();
        let samples: Vec<f64> = (0..100_000).map(|_| normal.sample(&mut rng)).collect();

        let welch = WelchBuilder::new(samples)
            .with_fs(4.)
            .with_dft_size(128)
            .with_overlap_size(1)
            .with_segment_size(8)
            .with_normalization(Normalization::Density)
            .build();
        assert_relative_eq!(
            welch.periodogram().sum::<f64>()
                / welch.periodogram().collect::<Vec<f64>>().len() as f64,
            0.21985243050737127,
            epsilon = 1e-2
        );
    }
}
