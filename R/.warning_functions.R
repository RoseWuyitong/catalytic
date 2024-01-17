.large_number_warning <- function(num, threshold, suggestion) {
  if (num > threshold) warning(suggestion)
}

.small_number_warning <- function(num, threshold, suggestion) {
  if (num < threshold) warning(paste(suggestion))
}

.not_same_stop <- function(num1, num2, suggestion) {
  if (num1 != num2) stop(suggestion)
}

.not_in_stop <- function(obj, expected_objs, suggestion) {
  if (obj %in% expected_objs) stop(suggestion)
}

.not_single_str_stop <- function(obj, suggestion) {
  if (!is.character(obj) | length(obj) != 1) stop(suggestion)
}
