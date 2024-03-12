posterior_predict.cat_glm <- function(cat_glm,
                              newdata,
                              ...) {

  return(matrix_multiply(standardize_data(newdata), cat_glm$est_beta))
}
