#' Discrete Fourier Transform
#'
#' @param data vector containing a numeric time series in the time-intensity
#' domain to be transformed to the frequency-amplitude domain
#' @param freq_seq vector of
#'
#' @return
#' @export
#'
#' @examples
dft <- function(data, freq_seq = NULL, sampling_rate = NULL) {
  n <- length(data)
  if(is.null(freq_seq)) {
    freq_seq <- 1:n
  }
  else if(!is.null(sampling_rate)) {
    freq_seq <- freq_seq * (n/sampling_rate)
  }
  index <- 1:length(data)
  per_freq <- lapply(freq_seq, function(freq) {
      angles <- freq * index / n * 2 * pi
      return(get_amplitude(data*cos(angles), data*sin(angles)))
    })

  if(!is.null(sampling_rate)) {
    freq_seq = freq_seq / (n / sampling_rate)
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
get_peaks <- function(data, decline=0.3, window=15) {
  shape <- diff(sign(diff(data, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - window + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + window + 1
    w <- ifelse(w < length(data), w, length(data))
    if(all(data[c(z:i, (i + 2):w)] <= (data[i + 1]*(1-decline)))) {
      return(i + 1)
      }
    else {
      return(numeric(0))
    }
  })
  pks <- unlist(pks)
  return(pks)
}

#' Title
#'
#' @param panel
#'
#' @return
#'
#' @examples
fft <- function(panel) {

}

#' Title
#'
#' @param panel
#'
#' @return
#'
#' @examples
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
  return(2*sqrt(mean(x)^2 +  mean(y)^2))
}
