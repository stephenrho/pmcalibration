#' fits a calibration curve via glm
#'
#' @param y binary outcome
#' @param p predicted probabilities
#' @param X design matrix created via \code{reg_spline_X} (rows of \code{X} should align with elements of \code{p})
#' @param Xp design matrix for to-be-plotted probability
#' @param save_data whether to save the data elements in the returned object
#' @param save_mod whether to save the model in the returned object
#'
#' @returns list of class \code{glm_cal}
#' @keywords internal
#' @export
glm_cal <- function(y, p, X, Xp, save_data = T, save_mod = T){

  # fit the calibration curve model
  d <- data.frame(y, X)

  mod <- glm(y ~ ., data = d, family = binomial(link="logit"))

  p_c <- predict(mod, type = "response")

  if (!is.null(Xp)){
    p_c_plot <- predict(mod, newdata = Xp, type = "response")
  } else{
    p_c_plot <- NULL
  }

  out <- list(
    y = if (save_data) y else NULL,
    p = if (save_data) p else NULL,
    X = if (save_data) X else NULL,
    Xp = if (save_data) Xp else NULL,
    p_c = p_c, # predict(mod, type = "response"),
    metrics = cal_metrics(p, p_c),
    p_c_plot = p_c_plot, # predict(mod, newdata = Xp, type = "response"),
    model = if (save_mod) mod else NULL
  )

  class(out) <- "glm_cal"
  return(out)
}
