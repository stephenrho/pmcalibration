#' calibration curve via \code{lowess}
#'
#' @description
#' uses \code{lowess} with \code{iter} = 0; as done in the \code{rms} package
#'
#' @param y binary outcome
#' @param p predicted probabilities
#' @param x predictor (could be transformation of \code{p})
#' @param xp values for plotting (same scale as \code{x})
#' @param save_data whether to save y, p, x, xp in the returned object
#'
#' @returns list of class \code{lowess_cal}
#' @keywords internal
#' @export
lowess_cal <- function(p, y, x, xp, save_data = T){

  fit <- lowess(x = x, y = y, iter = 0)

  p_c <- predict_lowess(fit = fit, x = x)
  if (!is.null(xp)){
    p_c_plot <- predict_lowess(fit = fit, x = xp)
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
    p_c_plot = p_c_plot
  )

  class(out) <- "lowess_cal"
  return(out)
}

#' @keywords internal
#' @export
predict_lowess <- function(fit, x){
  # adapted from rms:::calibrate.default
  pred = approx(x = fit$x, y = fit$y, xout=x, ties=function(xx)xx[1], rule=2)$y
  # rule = 2 means that if x is outside range(fit$x) the closest value is returned
  # see ?approx
  # this will happen on boostrap resamples for xp
  # better to return NAs? not sure...
  return(pred)
}

