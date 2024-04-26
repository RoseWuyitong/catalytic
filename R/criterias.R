# Based on paper https://arxiv.org/pdf/1307.5928.pdf


predict_error <- function(X,
                          Y,
                          family,
                          estimate_method = c(
                            "logarithmic_error",
                            "square_error",
                            "classification_error"
                          ),
                          estimate_coefficients = NULL,
                          estimate_Y = NULL,
                          estimate_cat_glm_model = NULL) {
  fam <- get_family(family)
  return()
}
