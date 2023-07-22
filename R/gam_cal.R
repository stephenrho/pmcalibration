#' fits a calibration curve via gam
#'
#' @param y binary outcome
#' @param p predicted probabilities
#' @param x predictor (could be transformation of \code{p})
#' @param xp values for plotting (same scale as \code{x})
#' @param save_data whether to save the data elements in the returned object
#' @param save_mod whether to save the model in the returned object
#' @param ... additional arguments for \code{mgcv::gam} and \code{mgcv::s}
#'
#' @returns list of class \code{gam_cal}
#' @keywords internal
#' @export
gam_cal <- function(y, p, x, xp, save_data = T, save_mod = T, ...){

  dots <- list(...)
  if ("bs" %in% names(dots)) bs <- dots[['bs']] else bs <- "tp"
  if ("k" %in% names(dots)) k <- dots[['k']] else k <- -1
  if ("fx" %in% names(dots)) fx <- dots[['fx']] else fx <- FALSE
  if ("method" %in% names(dots)) method <- dots[['method']] else method <- "GCV.Cp"

  # fit the calibration curve model
  d <- data.frame(y, x)

  mod <- mgcv::gam(y ~ s(x, k = k, fx = fx, bs = bs), data = d,
                   family = binomial(link="logit"),
                   method = method)

  p_c <- as.vector(predict(mod, type = "response"))

  if (!is.null(xp)){
    p_c_plot <- as.vector(predict(mod, newdata = data.frame(x = xp), type = "response"))
  } else{
    p_c_plot <- NULL
  }

  out <- list(
    y = if (save_data) y else NULL,
    p = if (save_data) p else NULL,
    x = if (save_data) x else NULL,
    xp = if (save_data) xp else NULL,
    p_c = p_c, # predict(mod, type = "response"),
    metrics = cal_metrics(p, p_c),
    p_c_plot = p_c_plot, # predict(mod, newdata = Xp, type = "response"),
    model = if (save_mod) mod else NULL,
    smooth_args = list(
      smooth = "gam",
      bs = bs,
      k = k,
      fx = fx,
      method = method
    )
  )

  class(out) <- "gam_cal"
  return(out)
}

