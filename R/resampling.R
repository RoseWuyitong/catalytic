#' resampling data
#'
#' @param original_data numeric or complex or logical vectors.
#' Data used for resampling.
#' @param resample_data_size integer.
#' Number of resampled output data.
#' @param resample_method str.
#' Choose from (coordinate,
#' symmetrizing,
#' smoothing,
#' flattering,
#' simple).
#' @param smoothing_noise_variance numeric.
#' The variance used for smoothing resample method only.
#'
#' @return matrix.
#' Resampled data.
#'
#' @export
#'
#' @examples
#' ## Pass a vector to resampling, and it will transform it into matrix before further process
#' resampling(1:5, resample_method = "coordinate")
#'
#' ## Pass the same input, and it will output the memerised result
#' resampling(
#'   1:5,
#'   resample_method = "coordinate"
#' ) == resampling(
#'   1:5,
#'   resample_method = "coordinate"
#' )
#'
#' ## Pass the different input, and it will output the processed result
#' resampling(1:5,
#'   resample_data_size = 5,
#'   resample_method = "coordinate"
#' )
resampling <- memoise::memoise(function(original_data,
                                        resample_data_size = 10,
                                        resample_method = "simple",
                                        smoothing_noise_variance = 1) {
  original_data <- as.matrix(original_data)
  resampled_data <- switch(resample_method,
    "coordinate" = apply(original_data,
      2,
      sample,
      size = resample_data_size,
      replace = TRUE
    ),
    "symmetrizing" = apply(
      original_data,
      2,
      function(x) {
        sample(c(x, 2 * mean(original_data) - x),
          resample_data_size,
          replace = TRUE
        )
      }
    ),
    "smoothing" = original_data + stats::rnorm(
      length(original_data),
      0,
      sqrt(smoothing_noise_variance)
    ),
    "flattering" = apply(
      original_data,
      2,
      function(x) {
        sample(c(x, sample(c(unique(original_data)), 1)),
          resample_data_size,
          replace = TRUE
        )
      }
    ),
    "simple" = sample(original_data,
      resample_data_size,
      replace = TRUE
    ),
    NA
  )

  return(resampled_data)
})
