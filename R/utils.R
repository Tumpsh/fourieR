#' Discrete Fourier Transform
#'
#' @param data vector containing a numeric time series in the time-intensity
#' domain to be transformed to the frequency-amplitude domain
#' @param freq_seq vector of frequencies to include in our transformed data. If
#' null (default), will be 1:length(data). If sampling_rate is not null,
#' provided frequencies will be consider as per sampling_rate samples.
#' @param sampling_rate single number, indicating number of samples per unit.
#' Unit should be the same as the denominator for the units of frequency being
#' considered.
#'
#' @return data.frame containing frequencies and their amplitudes in two
#' columns.
#'
#' @export
#'
#' @examples
#' audio_file <- system.file("extdata/AudioSampleComplex.mp3",
#' package = "fourieR")
#' library(tuneR)
#' audio_wave <- tuneR::readMP3(audio_file)
#' wave_form <- audio_wave@@left
#' wave_freq <- dft(wave_form, freq_seq = 1:2000,
#'                  sampling_rate = audio_wave@@samp.rate)
#' peak_freqs <- get_peaks(wave_freq$amp)
#' cleaned_wave_form <- inverse_dft(wave_freq[peak_freqs, ], 3,
#'                                  audio_wave@@samp.rate)
#' cleaned_wave <- Wave(cleaned_wave_form, cleaned_wave_form,
#'                     audio_wave@@samp.rate)
#' ## writeWave(cleaned_wave, "cleaned_wav.wav")
#' ## writeWave(audio_wave, "original_wav.wav")
dft <- function(data, freq_seq = NULL, sampling_rate = NULL) {
  n <- length(data)
  if (is.null(freq_seq)) {
    freq_seq <- 1:n
  }
  else if (!is.null(sampling_rate)) {
    freq_seq <- freq_seq * (n / sampling_rate)
  }
  index <- 1:length(data)
  per_freq <- lapply(freq_seq, function(freq) {
    angles <- freq * index / n * 2 * pi
    return(get_amplitude(data * cos(angles), data * sin(angles)))
  })

  if (!is.null(sampling_rate)) {
    freq_seq <- freq_seq / (n / sampling_rate)
  }
  return(data.frame(freq = freq_seq, amp = unlist(per_freq)))
}

#' Get Peaks
#'
#' Implements a local maximum finding algorithm, with certain tolerances for
#' how sharp the peak is, and how spread apart peaks must be.
#'
#' @param data vector of numeric data
#' @param decline optional, between (0, 1]. Fraction of point P which nearby
#' points are expected be lower than if P is a true peak (e.g. If the
#' adjacent points to P are lower than 0.8*P, then we will consider P a peak)
#' @param window optional, number of points on either side surrounding a point P
#' which are expected to be lower than decline*P if P is a true peak.
#'
#' @return Return the indices within data of the local maxima discovered
#' @export
#'
#' @examples
#' data <- c(1:5, 20, 6:1, rep(5, 5), 6, rep(5,5))
#' ## Should only return 6, since the element, 20, at the 6th index, has at
#' ## least 4 elements on either side which are at least 20% lower. Even though
#' ## element 18 also is a local maximum, it is not sufficiently higher than its
#' ## neighbors
#' get_peaks(data, 0.2, 4)
get_peaks <- function(data, decline = 0.3, window = 15) {
  shape <- diff(sign(diff(data, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i) {
    z <- i - window + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + window + 1
    w <- ifelse(w < length(data), w, length(data))
    if (all(data[c(z:i, (i + 2):w)] <= (data[i + 1] * (1 - decline)))) {
      return(i + 1)
    }
    else {
      return(numeric(0))
    }
  })
  pks <- unlist(pks)
  return(pks)
}

#' Inverse Discrete Fourier Transform
#'
#' @param freq_data data.frame containing columns "freq" and "amp" with
#' frequencies and their amplitudes for the inverse transform.
#' @param len length of the output waveform, in the unit of interest, or simply
#' number of data points, if no sampling_rate is given.
#' @inheritParams dft
#'
#' @return vector containing the intensity-time domain time-series.
#' @export
#'
#' @examples
#' audio_file <- system.file("extdata/AudioSampleComplex.mp3",
#' package = "fourieR")
#' library(tuneR)
#' audio_wave <- tuneR::readMP3(audio_file)
#' wave_form <- audio_wave@@left
#' wave_freq <- dft(wave_form, freq_seq = 1:2000,
#'                  sampling_rate = audio_wave@@samp.rate)
#' peak_freqs <- get_peaks(wave_freq$amp)
#' cleaned_wave_form <- inverse_dft(wave_freq[peak_freqs, ], 3,
#'                                  audio_wave@@samp.rate)
#' cleaned_wave <- Wave(cleaned_wave_form, cleaned_wave_form,
#'                     audio_wave@@samp.rate)
#' ## writeWave(cleaned_wave, "cleaned_wav.wav")
#' ## writeWave(audio_wave, "original_wav.wav")
inverse_dft <- function(freq_data, len, sampling_rate = NULL) {
  samples <- 1:round(len * ifelse(is.null(sampling_rate), 1, sampling_rate))
  res <- rep(0, length(samples))
  if (is.null(sampling_rate)) {
    sampling_rate <- len
  }
  for (row_num in 1:nrow(freq_data)) {
    pair <- freq_data[row_num, ]
    res <- res + pair[, "amp"] * sin(2 * pi * samples * pair[, "freq"] /
      (sampling_rate))
  }
  return(res)
}

fft <- function(panel) {

}

nonunif_fft <- function(panel) {

}

#' Get Amplitude of the Center of Mass of the Wrapped Signal
#'
#' For a certain signal at a frequency, the signal is "wrapped" around the
#' origin at a specific frequency. The magnitude of the center of mass of that
#' wrapped function is half of the amplitude of that frequency in the transform.
#'
#' @param x,y Vectors of X and Y coordinates
#'
#' @return amplitude of the relevant frequency used to produce the X and Y
#' coordinates in the original signal.
#'
get_amplitude <- function(x, y) {
  return(2 * sqrt(mean(x)^2 + mean(y)^2))
}
