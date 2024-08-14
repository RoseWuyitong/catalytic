#' Title
#'
#' @param cat_glm_
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
predict.cat_glm <- function(cat_glm_, ...) {
  stats::predict.glm(cat_glm_$model,
    type = if (cat_glm_$cat_init$family == "binomial") "response" else "link",
    ...
  )
}


#' Title
#'
#' @param bayes_
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
predict.bayes <- function(bayes_,
                          newdata = NULL) {
  get_glm_mean(
    bayes_$cat_init$family,
    if (is.null(newdata)) bayes_$cat_init$obs_x else newdata[, colnames(bayes_$cat_init$obs_x)],
    coef(bayes_)
  )
}
