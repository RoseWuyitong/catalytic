#' generate data from distribution
#'
#' @param from_dist str.
#' The distribution for generating data, choose from
#' (bernoulli,
#' gaussian,
#' uniform)
#' @param sample_size integer.
#' The number of output sample data
#' @param dist_mean double.
#' The mean of the distribution for generating data
#' @param dist_sd double.
#' The standard deviation of the distribution for generating data
#'
#' @return matrix. sample data
#' @export
#'
#' @examples
#' ## Generate 10 data from Gaussian distribution with mean 0 and sd 1
#' generate_data_from_dist(from_dist = "gaussian")
#'
#' ## Pass the same input, and it will output the memerised result
#' generate_data_from_dist(
#'   from_dist = "gaussian"
#' ) == generate_data_from_dist(
#'   from_dist = "gaussian"
#' )
#'
#' ## Pass the different input, and it will output the processed result
#'
#' generate_data_from_dist(
#'   from_dist = "gaussian",
#'   sample_size = 5
#' )
generate_data_from_dist <- memoise::memoise(function(from_dist = "uniform",
                                                     sample_size = 10,
                                                     dist_mean = 0,
                                                     dist_sd = 1) {
  result <- switch(from_dist,
    "bernoulli" = stats::rbinom(
      sample_size,
      1,
      meadist_meann_
    ),
    "gaussian" = stats::rnorm(
      sample_size,
      dist_mean,
      dist_sd
    ),
    "uniform" = rep(dist_mean, sample_size)
  )

  return(as.matrix(result))
})
