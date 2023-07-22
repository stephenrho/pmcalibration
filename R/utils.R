#' Summarize a pmcalibration object
#'
#' @param x object created with \code{pmcalibration}
#' @param conf_level width of the confidence interval (0.95 gives 95\% CI). Ignored if call to \code{pmcalibration} didn't request confidence intervals
#'
#' @return prints a summary of calibration metrics. Returns a list of two tables: \code{metrics} and \code{plot}
#' @export
summary.pmcalibration <- function(x, conf_level = .95){

  checkmate::assert_double(conf_level, len = 1)

  probs <- c((1 - conf_level)/2, 1 - (1 - conf_level)/2)

  m_tab <- data.frame(Estimate = x$metrics)

  if (!is.null(x$metrics.samples)){
    m_ci <- t(apply(x$metrics.samples, MARGIN = 2, FUN = quantile, probs = probs, na.rm = T))
    colnames(m_ci) <- c("lower", "upper")
    m_tab <- cbind(m_tab, m_ci)
  }

  plot_tab <- data.frame(p = x$plot$p, p_c = x$plot$p_c_plot)

  if (!is.null(x$plot$plot.samples)){
    p_ci <- t(apply(x$plot$plot.samples, MARGIN = 2, FUN = quantile, probs = probs, na.rm = T))
    colnames(p_ci) <- c("lower", "upper")
    plot_tab <- cbind(plot_tab, p_ci)
  }

  #print(m_tab, digits = 3)

  out <- list(
    metrics = m_tab,
    plot = plot_tab,
    smooth = x$smooth,
    ci = x$ci,
    conf_level = conf_level,
    n = x$n,
    smooth_args = x$smooth_args,
    logitp = x$logitp
    )

  class(out) <- "pmcalibratesummary"

  return(out)
}

#' @export
print.pmcalibratesummary <- function(x, digits = 2){
  smoothtext = list("none" = "no smooth",
                    "rcs" = "a restricted cubic spline (see ?rms::rcs)",
                    "ns" = "a natural cubic spline (see ?splines::ns)",
                    "bs" = "a B-spline (see ?splines::bs)",
                    "gam" = "a generalized additive model (see ?mgcv::s)"
                    )

  citext <- list("boot" = "bootstrap resampling",
                 "sim" = "simulation based inference"
                 )

  transp <- if (x$logitp) "logit transformed" else "untransformed"

  cat("Calibration metrics based on a calibration curve estimated via",
      smoothtext[[x$smooth]], "using", transp, "predicted probabilities.\n\n")

  print(x$metrics, digits = digits)

  if (x$ci %in% names(citext)){
    cat("\n")
    cat(sprintf("%.0f%%", x$conf_level*100),
        "confidence intervals calculated via",
        citext[[x$ci]], "with", x$n, "replicates.\n")
  }

  invisible(x)
}

#' @export
print.pmcalibration <- function(x, digits = 2, conf_level = .95) {
  print(summary(x, conf_level = conf_level), digits = digits)
}


#' Extract plot data from \code{pmcalibration} object
#'
#' @param x \code{pmcalibration} object
#' @param conf_level width of the confidence interval (0.95 gives 95\% CI). Ignored if call to \code{pmcalibration} didn't request confidence intervals
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
get_cc <- function(x, conf_level = .95){
  # cc <- data.frame(
  #   p = x$plot$p,
  #   p_c = x$plot$p_c_plot,
  #   lower = x$plot$conf.int[, 1],
  #   upper = x$plot$conf.int[, 2]
  #   )

  cc <- summary(x, conf_level = conf_level)$plot
  return(cc)
}

#' Plot a calibration curve
#'
#' @description
#' This is for a quick and dirty calibration curve plot.
#' Alternatively you can use \code{get_cc()} to get the data required to plot the calibration curve.
#'
#' @param x a \code{pmcalibration} calibration curve
#' @param conf_level width of the confidence interval (0.95 gives 95\% CI). Ignored if call to \code{pmcalibration} didn't request confidence intervals
#' @param ... other args for \code{plot()} (\code{lim} and \code{lab} can be specified)
#'
#' @export
plot.pmcalibration <- function(x, conf_level = .95, ...){

  dots <- list(...)

  # pdat <- summary.pmcalibration(x, conf_level = conf_level)
  # pdat <- pdat$plot
  pdat <- get_cc(x, conf_level = conf_level)

  if ("xlim" %in% names(dots)) xlim <- dots[['xlim']] else xlim <- c(0,1)
  if ("ylim" %in% names(dots)) ylim <- dots[['ylim']] else ylim <- c(0,1)
  if ("xlab" %in% names(dots)) xlab <- dots[['xlab']] else xlab <- "Predicted Probability"
  if ("ylab" %in% names(dots)) ylab <- dots[['ylab']] else ylab <- "Estimated Probability"

  plot(x = pdat$p, y = pdat$p_c, type="l",
       xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

  if ("lower" %in% colnames(pdat)){
    lines(x = pdat$p, y = pdat$lower, lty=2)
    lines(x = pdat$p, y = pdat$upper, lty=2)
  }
}


#' @keywords internal
#' @export
logit <- binomial()$linkfun

#' @keywords internal
#' @export
invlogit <- binomial()$linkinv
