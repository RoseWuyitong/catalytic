#' generate synthetic data
#'
#' @param observed_covariate numeric or complex or logical vectors.
#' Observed data, x value.
#' @param observed_response numeric or complex or logical vectors.
#' Observed response, y value.
#' @param covariate_generate_from str.
#' Choose from (resampling, distribution), resampling is resample the original data,
#' while distribution is generate data from distribution instead of original data.
#' @param response_generatie_from str.
#' Choose from (resampling, distribution), resampling is resample the original data,
#' while distribution is generate data from distribution instead of original data.
#' @param synthetic_data_size integer.
#' The number of synthetic data.
#' @param covariate_resample_method str.
#' Choose from (coordinate,
#' symmetrizing,
#' smoothing,
#' flattering,
#' simple).
#' @param response_resample_method str.
#' Choose from (coordinate,
#' symmetrizing,
#' smoothing,
#' flattering,
#' simple).
#' @param covariate_smoothing_noise_variance numeric.
#' The variance used for smoothing resample method only.
#' @param response_smoothing_noise_variance numeric.
#' The variance used for smoothing resample method only.
#' @param covariate_distribution str.
#' The distribution for generating data, choose from
#' (bernoulli,
#' gaussian,
#' uniform)
#' @param response_distribution str.
#' The distribution for generating data, choose from
#' (bernoulli,
#' gaussian,
#' uniform)
#' @param covariate_distribution_mean double.
#' The mean of the distribution for generating data
#' @param response_distribution_mean double.
#' The mean of the distribution for generating data
#' @param covariate_distribution_sd double.
#' The standard deviation of the distribution for generating data
#' @param response_distribution_sd double.
#' The standard deviation of the distribution for generating data
#'
#' @return rlang object.
#' Contains observed data and synthetic data
#' env(obs_x, obs_y, syn_x, syn_y).
#' @export
#'
#' @examples
#' ## Generate synthetic covariate with resample covariate with coordinate, and synthetic response with Gaussian distribution
#' data <- synthetic_data_generation(
#' matrix(1:9, nrow = 3, ncol = 3), c(10, 20, 30),
#' covariate_resample_method="coordinate",
#' response_distribution="gaussian")
#'
#' ## Observed data
#' data$obs_x
#' data$obs_y
#'
#' ## Synthetic data
#' data$syn_x
#' data$syn_y
#'
#' ## Pass the same input, and it will output the memerised result
#' synthetic_data_generation(
#' matrix(1:9, nrow = 3, ncol = 3), c(10, 20, 30),
#' covariate_resample_method="coordinate",
#' response_distribution="gaussian") == synthetic_data_generation(
#' matrix(1:9, nrow = 3, ncol = 3), c(10, 20, 30),
#' covariate_resample_method="coordinate",
#' response_distribution="gaussian")
#'
#' ## Pass the different input, and it will output the processed result
#' data_ <- synthetic_data_generation(
#' matrix(1:9, nrow = 3, ncol = 3), c(10, 20, 30),
#' covariate_resample_method="coordinate",
#' response_distribution="gaussian",
#' synthetic_data_size=5)
#' data_$syn_x
#' data_$syn_y
synthetic_data_generation <- memoise::memoise(function(observed_covariate,
                                                       observed_response,
                                                       covariate_generate_from = "resampling",
                                                       response_generate_from = "distribution",
                                                       synthetic_data_size = 10,
                                                       covariate_resample_method = "simple",
                                                       response_resample_method = "simple",
                                                       covariate_smoothing_noise_variance = 1,
                                                       response_smoothing_noise_variance = 1,
                                                       covariate_distribution = "uniform",
                                                       response_distribution = "uniform",
                                                       covariate_distribution_mean = 0,
                                                       response_distribution_mean = 0,
                                                       covariate_distribution_sd = 1,
                                                       response_distribution_sd = 1) {
  data <- rlang::env(
    obs_x = as.matrix(observed_covariate),
    obs_y = as.matrix(observed_response)
  )


  data$syn_x <- switch(covariate_generate_from,
    "resampling" = resampling(
      original_data = observed_covariate,
      resample_data_size = synthetic_data_size,
      resample_method = covariate_resample_method,
      smoothing_noise_variance = covariate_smoothing_noise_variance
    ),
    "distribution" = generate_data_from_dist(
      from_dist = covariate_distribution,
      sample_size = synthetic_data_size,
      dist_mean = covariate_distribution_mean,
      dist_sd = covariate_distribution_sd
    )
  )

  data$syn_y <- switch(response_generate_from,
    "resampling" = resampling(
      original_data = observed_response,
      resample_data_size = synthetic_data_size,
      resample_method = response_resample_method,
      smoothing_noise_variance = response_smoothing_noise_variance
    ),
    "distribution" = generate_data_from_dist(
      from_dist = response_distribution,
      sample_size = synthetic_data_size,
      dist_mean = response_distribution_mean,
      dist_sd = response_distribution_sd
    )
  )

  return(data)
})
