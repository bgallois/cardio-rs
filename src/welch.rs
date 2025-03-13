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
#![cfg(feature = "std")]

use num::{Complex, Float};
use rustfft::{Fft, FftDirection, FftNum};
use std::{cmp, f64::consts::PI};

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
        T: Float + Copy + core::fmt::Debug + FftNum + std::iter::Sum,
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

/// Represents the result of Welch's method for power spectral density estimation.
///
/// Welch's method reduces noise in power spectral estimates by averaging
/// multiple periodograms computed from overlapping segments.
///
/// # Type Parameters
/// * `T` - A numeric type that supports floating-point operations.
pub struct Welch<T> {
    pub periodogram: Vec<T>,
    pub frequencies: Vec<T>,
}

/// Builder for constructing a Welch power spectral density estimator.
///
/// This struct provides a flexible way to configure the parameters for Welch's method,
/// such as segment size, overlap, and FFT size.
///
/// # Type Parameters
/// * `T` - A floating-point type that supports FFT operations.
pub struct WelchBuilder<T> {
    segment_size: usize,
    dft_size: usize,
    overlap_size: usize,
    fs: T,
    signal: Vec<T>,
}

impl<T: Float + Copy + core::fmt::Debug + FftNum + std::iter::Sum + std::ops::AddAssign>
    WelchBuilder<T>
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
            let chunk: Vec<T> =
                self.signal[i..cmp::min(i + self.segment_size, self.signal.len())].to_vec();
            let chunk_mean = chunk.iter().copied().sum::<T>() / T::from(chunk.len()).unwrap();
            let chunk = chunk
                .iter()
                .zip(window.weights.iter())
                .map(|(&a, &b)| (a - chunk_mean) * b)
                .collect::<Vec<T>>();

            let fft = rustfft::algorithm::Radix4::new(self.dft_size, FftDirection::Forward);

            let mut buffer = vec![
                Complex {
                    re: T::from(0).unwrap(),
                    im: T::from(0).unwrap(),
                };
                self.dft_size
            ];
            chunk.iter().enumerate().for_each(|(i, j)| {
                buffer[i].re = *j;
            });
            fft.process(&mut buffer);

            let pdg: Vec<T> = buffer
                .iter()
                .take(self.dft_size / 2)
                .map(|i| i.norm_sqr())
                .collect();
            periodogram
                .iter_mut()
                .zip(pdg.iter())
                .for_each(|(a, b)| *a += *b);

            i += self.segment_size - self.overlap_size;
            n_segments += 1;
        }
        let frequencies: Vec<T> = (0..self.dft_size / 2)
            .map(|i| T::from(i).unwrap() * self.fs / T::from(self.dft_size).unwrap())
            .collect();
        let window_power: T = window.weights.iter().map(|&w| w * w).sum();
        let norma = (window_power * T::from(n_segments).unwrap() * self.fs).recip();
        let periodogram = periodogram.iter().map(|&p| p * norma).collect();

        Welch {
            periodogram,
            frequencies,
        }
    }
}
