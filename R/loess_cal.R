#' calibration curve via \code{loess}
#'
#' @param y binary outcome
#' @param p predicted probabilities
#' @param x predictor (could be transformation of \code{p})
#' @param xp values for plotting (same scale as \code{x})
#' @param save_data whether to save y, p, x, xp in the returned object
#' @param save_mod whether to save the model in the returned object
#'
#' @returns list of class \code{loess_cal}
#' @keywords internal
#' @export
loess_cal <- function(p, y, x, xp, save_data = T, save_mod = T){

  mod <- loess(y ~ x)
  # TODO
  # explore loess options more fully...
  # control = loess.control(surface="direct")
  # would allow extrapolation on bs resamples but seems to take longer...

  p_c <- predict(mod, newdata = x)
  if (!is.null(xp)){
    p_c_plot <- predict(mod, newdata = xp)
  } else{
    p_c_plot <- NULL
  }

  out <- list(
    y = if (save_data) y else NULL,
    p = if (save_data) p else NULL,
    x = if (save_data) x else NULL,
    xp = if (save_data) xp else NULL,
    p_c = p_c,
    metrics = cal_metrics(p, p_c),
    p_c_plot = p_c_plot,
    model = if (save_mod) mod else NULL,
    smooth_args = NULL
  )

  class(out) <- "loess_cal"
  return(out)
}
