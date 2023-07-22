#' fits a calibration curve via glm
#'
#' @param y binary outcome
#' @param p predicted probabilities
#' @param x predictor (could be transformation of \code{p})
#' @param xp values for plotting (same scale as \code{x})
#' @param smooth 'rcs', 'ns', 'bs', or 'none'
#' @param save_data whether to save the data elements in the returned object
#' @param save_mod whether to save the model in the returned object
#' @param pw save pointwise standard errors for plotting
#'
#' @returns list of class \code{glm_cal}
#' @keywords internal
#' @export
glm_cal <- function(y, p, x, xp, smooth, save_data = T, save_mod = T, pw = F, ...){

  XX <- reg_spline_X(x = x, xp = xp, smooth = smooth, ...)
  X <- XX$X
  Xp <- XX$Xp

  # fit the calibration curve model
  d <- data.frame(y, X)

  mod <- glm(y ~ ., data = d, family = binomial(link="logit"))

  p_c <- predict(mod, type = "response")

  if (!is.null(Xp)){
    if (pw){
      p_c_p <- predict(mod, newdata = Xp, type = "response", se.fit = T)
      p_c_plot <- p_c_p$fit
      p_c_plot_se <- p_c_p$se
    } else{
      p_c_plot <- predict(mod, newdata = Xp, type = "response")
      p_c_plot_se <- NULL
    }
  } else{
    p_c_plot <- NULL
    p_c_plot_se <- NULL
  }

  out <- list(
    y = if (save_data) y else NULL,
    p = if (save_data) p else NULL,
    x = if (save_data) x else NULL,
    xp = if (save_data) xp else NULL,
    X = if (save_data) X else NULL,
    Xp = if (save_data) Xp else NULL,
    p_c = p_c,
    metrics = cal_metrics(p, p_c),
    p_c_plot = p_c_plot,
    p_c_plot_se = p_c_plot_se,
    model = if (save_mod) mod else NULL,
    smooth_args = XX$smooth_args
  )

  class(out) <- "glm_cal"
  return(out)
}
