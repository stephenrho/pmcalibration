#' Check if method and smooth are allowed together (not used anymore)
#'
#' @keywords internal
#' @export
check_method_smooth <- function(method, smooth){

  allowed <- list("glm" = c("none", "ns", "bs", "rcs"),
                  "gam" = c("s"))

  if (method %in% c("glm", "gam")){
    if (!smooth %in% allowed[[method]]){
      stop(sprintf("smooth = %s is not allowed with method = %s", smooth, method))
    }
  } else if (method %in% c("loess", "lowess") & !is.null(smooth)){
    warning("smooth is ignored for loess and lowess")
  }
}

#' Get confidence interval from objects created with \code{boot()} or \code{simb()}
#'
#' @keywords internal
#' @export
get_ci <- function(b, conf_level = .95){

  probs <- c((1 - conf_level)/2, 1 - (1 - conf_level)/2)

  m <- do.call(rbind, lapply(b, function(x) x$metrics))

  if (any(is.na(m))){
    warning("Bootstrap resamples contain NAs. This is probably due to a loess not extrapolating for to-be-plotted values...")
  }

  m_ci <- t(apply(m, MARGIN = 2, FUN = quantile, probs = probs, na.rm = T))

  if (!is.null(b[[1]]$p_c_plot)){
    pplot <- do.call(rbind, lapply(b, function(x) x$p_c_plot))
    if (any(is.na(pplot))){
      warning("Bootstrap resamples contain NAs. This is probably due to a loess not extrapolating for to-be-plotted values...")
    }

    pplot_ci <- t(apply(pplot, MARGIN = 2, FUN = quantile, probs = probs, na.rm = T))
  } else{
    pplot_ci <- NULL
  }

  cis <- list(metrics = m_ci, p_c_plot = pplot_ci)

  return(cis)
}

#' Plot a calibration curve
#'
#' @description
#' This is for a quick and dirty calibration curve plot.
#' Alternatively you can use \code{get_cc()} to get the data required to plot the calibration curve.
#'
#' @param x a \code{pmcalibration} calibration curve
#' @param ... other args for plot()
#'
#' @export
plot.pmcalibration <- function(x, ...){

  dots <- list(...)

  p <- x$plot$p
  p_c <- x$plot$p_c_plot
  ci <- x$plot$conf.int

  if ("xlim" %in% names(dots)) xlim <- dots[['xlim']] else xlim <- c(0,1)
  if ("ylim" %in% names(dots)) ylim <- dots[['ylim']] else ylim <- c(0,1)
  if ("xlab" %in% names(dots)) xlab <- dots[['xlab']] else xlab <- "Predicted Probability"
  if ("ylab" %in% names(dots)) ylab <- dots[['ylab']] else ylab <- "Estimated Probability"

  plot(x = p, y = p_c, type="l", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  lines(x = p, y = ci[,1], lty=2)
  lines(x = p, y = ci[,2], lty=2)

}

#' Extract plot data from \code{pmcalibration} object
#'
#' @param x \code{pmcalibration} object
#'
#' @return data frame for plotting with 4 columns
#' \itemize{
#' \item{\code{p} - values for the x-axis (predicted probabilities - note these are *not* from your data and are only used for plotting)}
#' \item{\code{p_c} - probability impied by the calibration curve given \code{p}}
#' \item{\code{lower} and \code{upper} - bounds of the confidence interval
#' (type of CI and width determined by original call to \code{pmcalibration})}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' get_cc(cal) |>
#'  ggplot(aes(x = p, y = p_c, ymin=lower, ymax=upper)) +
#'  geom_abline(intercept = 0, slope = 1, lty=2) +
#'  geom_line() +
#'  geom_ribbon(alpha = 1/4)
#'  }
get_cc <- function(x){
  cc <- data.frame(
    p = x$plot$p,
    p_c = x$plot$p_c_plot,
    lower = x$plot$conf.int[, 1],
    upper = x$plot$conf.int[, 2]
    )
  return(cc)
}

#' @keywords internal
#' @export
logit <- binomial()$linkfun

#' @keywords internal
#' @export
invlogit <- binomial()$linkinv
