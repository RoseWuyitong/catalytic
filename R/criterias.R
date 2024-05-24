AIC <- function(...) {
  stats::AIC(...)
}

BIC <- function(...) {
  stats::BIC(...)
}

# anova
lrtest <- function(...) {
  lmtest::lrtest(randint1, randslope1, randslope2)
}

ICC <- function(...) {
  performance::icc(...)
}
