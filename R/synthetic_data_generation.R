data_env <- rlang::env(hash=TRUE)

synthetic_data_generation <- function(x
                                      )
{
return(x)
}


#' resampling data
#'
#' @param original_data data used for resampling
#' @param resample_data_size number of resampled output data
#' @param resample_method choose from coordinate_resampling,
#' symetrizing_resampling,
#' smoothing,
#' flattering_resampling,
#' simple_resampling
#' @param smoothing_noise_variance the variance used for
#' only smoothing resample method
#'
#' @return resampled data
#' @export
#'
#' @examples
#' resampling(c(1,2,3,4,5))
resampling <- function(original_data,
                       resample_data_size=10,
                       resample_method="coordinate_resampling",
                       smoothing_noise_variance=1
) {

  input_key <- as.character(rlang::hash(input_data))

  original_data <- as.data.frame(original_data)
  switch(resample_method,
         "coordinate_resampling"=apply(original_data,
                                       2,
                                       sample,
                                       size=resample_data_size,
                                       replace=TRUE),
         "symetrizing_resampling"=apply(original_data,
                                        2,
                                        function(x) sample(c(x, 2*mean(original_data)-x),
                                                           resample_data_size,
                                                           replace=TRUE)),
         "smoothing"=original_data+stats::rnorm(length(original_data),
                                         0,
                                         sqrt(smoothing_noise_variance)),
         "flattering_resampling"=apply(original_data,
                                       2,
                                       function(x) sample(c(x, sample(c(unique(original_data)), 1)),
                                                          resample_data_size)),
         "simple_resampling"=sample(original_data,
                                    resample_data_size,
                                    replace=TRUE),
         "NA")
}
