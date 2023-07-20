#' Create a calibration curve
#'
#' @description
#' Assess calibration of clinical prediction models (agreement between predicted and observed probabilities) via different smooths.
#'
#' @param y a binary outcome
#' @param p predicted probabilities from a clinical prediction model
#' @param method use regression splines (glm), generalized additive models (gam), or local regression smoothers (lowess, loess) to estimate calibration curve
#' @param smooth what smooth to use. For method = 'glm' and smooth = 'none' logistic calibration is performed
#' @param ci what kind of confidence intervals to compute? 'sim' = simulation based inference; \code{n} samples are taken from a multivariate normal distribution with mean vector = coef(mod) and variance covariance = vcov(model). 'boot' = bootstrap resampling with \code{n} replicates. \code{y} and \code{p} are sampled with replacement and calibration curve is reestimated. Calibration metrics are calculated using each simulation or boot sample. For both options percentile confidence intervals are returned.
#' @param conf_level width of the confidence interval (0.95 gives 95\% CI)
#' @param n number of simulations or bootstrap resamples
#' @param logitp fit calibration curve on logit transform of \code{p}
#' @param neval number of points (equally spaced between \code{min(p)} and \code{max(p)}) to evaluate for plotting (0 or NULL = no plotting)
#' @param ... additional arguments for particular smooths. For ci = 'boot' the user is able to run samples in parallel (using the \code{parallel} package) by specifying a \code{cores} argument
#'
#' @references Austin PC, Steyerberg EW. (2019) The Integrated Calibration Index (ICI) and related metrics for quantifying the calibration of logistic regression models. \emph{Statistics in Medicine}. 38, pp. 1–15. https://doi.org/10.1002/sim.8281
#' @references Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina M., Steyerberg E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. \emph{Journal of Clinical Epidemiology}, 74, pp. 167-176
#'
#' @returns a \code{pmcalibration} object containing calibration metrics and values for plotting
#' @export
pmcalibration <- function(y, p,
                          smooth=c("none", "ns", "bs", "rcs", "gam", "lowess", "loess"),
                          ci = c("sim", "boot", "none"), conf_level=.95,
                          n=1000, logitp = T, neval=100, ...){

  # TODO
  # - add ci option "pw" for pointwise (plot only)
  # - implement nveval = 0 | NULL --> dont save plot info (just set pplot to NULL?) (DONE)
  # - lowess_cal and loess_cal and associated boot methods (DONE)
  # - gam and associated methods (DONE)
  # - add print method for pmcalibration
  # - merge method and smooth - does it make sense to have separate? (DONE)
  # - function to extract plot data (get_cc()) (DONE)
  # - if 'knots' specified then use same knots on each boot sample (DONE except gam)
  # - allow Boundary.knots arg for ns/bs

  dots <- list(...)

  if (ci == "boot"){
    if ("cores" %in% names(dots)){
      cores <- dots[["cores"]]
    } else{
      cores <- 1
    }
  }

  chk::vld_compatible_lengths(y, p)

  # method <- match.arg(method)
  smooth <- match.arg(smooth)
  ci <- match.arg(ci)

  # check_method_smooth(method, smooth)

  if (!is.null(neval) && neval > 0){
    pplot <- seq(min(p), max(p), length.out = neval)
  } else{
    pplot <- NULL
  }

  if (logitp){
    xp <- if (!is.null(pplot)) logit(pplot) else NULL
    x <- logit(p)
  } else{
    x <- p; xp <- pplot
  }

  # fit cc
  if (smooth %in% c("none", "ns", "bs", "rcs")){
    # XX <- reg_spline_X(x = x, xp = xp, smooth = smooth, ...)
    # cal <- glm_cal(y = y, p = p, X = XX$X, Xp = XX$Xp, save_data = T, save_mod = T)
    cal <- glm_cal(y = y, p = p, x = x, xp = xp, smooth = smooth, save_data = T, save_mod = T, ...)
  } else if (smooth == "lowess"){
    cal <- lowess_cal(y = y, p = p, x = x, xp = xp, save_data = T)
  } else if (smooth == "loess"){
    cal <- loess_cal(y = y, p = p, x = x, xp = xp, save_data = T, save_mod = T)
  } else if (smooth == "gam"){
    cal <- gam_cal(y = y, p = p, x = x, xp = xp, save_data = T, save_mod = T, ...)
  }

  if (ci == "boot"){
    b.cal <- run_boots(cal, R = n, cores = cores)
    conf.int <- get_ci(b.cal)
  } else if (ci == "sim"){
    b.cal <- simb(cal, R = n)
    conf.int <- get_ci(b.cal)
  } else{
    conf.int <- NULL
  }

  out <- list(
    metrics = cal$metrics,
    conf.int = conf.int$metrics,
    plot = list(
      p = pplot,
      p_c_plot = cal$p_c_plot,
      conf.int = conf.int$p_c_plot
    ),
    smooth = smooth,
    ci = ci, conf_level = conf_level,
    n = n,
    smooth_ags = cal$smooth_args
  )

  class(out) <- "pmcalibration"

  return(out)

}
