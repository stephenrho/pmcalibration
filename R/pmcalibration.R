#' Create a calibration curve
#'
#' @description
#' Assess calibration of clinical prediction models (agreement between predicted and observed probabilities) via different smooths.
#'
#' @param y a binary outcome
#' @param p predicted probabilities from a clinical prediction model
#' @param smooth what smooth to use. Available options:
#' \itemize{
#' \item{'rcs' = restricted cubic spline using \code{rms::rcs}. Optional arguments for this smooth are \code{nk} (number of knots; defaults to 5) and \code{knots} (knot positions; set by \code{Hmisc::rcs.eval} if not specified)}
#' \item{'ns' = natural spline using \code{splines::ns}. Optional arguments are \code{df} (default = 6), \code{knots}, \code{Boundary.knots} (see \code{?splines::ns})}
#' \item{'bs' = B-spline using \code{splines::bs}. Optional arguments are \code{df} (default = 6), \code{knots}, \code{Boundary.knots} (see \code{?splines::bs})}
#' \item{'gam' = generalized additive model via \code{mgcv::gam} and \code{mgcv::s}. Optional arguments are \code{bs}, \code{k}, \code{fx}, \code{method} (see \code{?mgcv::gam} and  \code{?mgcv::s})}
#' \item{'lowess' = uses \code{lowess(x, y, iter = 0)} based on \code{rms::calibrate}}
#' \item{'loess' = uses \code{loess} with all defaults}
#' \item{'none' = logistic regression with single predictor variable (performs logistic calibration when \code{logitp = T})}
#' }
#' 'rcs', 'ns', 'bs', and 'none' are fit via \code{glm} and 'gam' is fit via \code{mgcv::gam} with \code{family = Binomial(link="logit")}
#' @param ci what kind of confidence intervals to compute?
#' \itemize{
#' \item{'sim' = simulation based inference; \code{n} samples are taken from a multivariate normal distribution with mean vector = coef(mod) and variance covariance = vcov(model).}
#' \item{'boot' = bootstrap resampling with \code{n} replicates. \code{y} and \code{p} are sampled with replacement and calibration curve is reestimated. If \code{knots} are specified the same knots are used for each resample (otherwise they are calculated using resampled x)}
#' \item{'pw' = pointwise confidence intervals calculated via the standard errors produced by relevant \code{predict} methods. Only for plotting curves - CIs are not produced for metrics (not available for smooth = 'lowess')}
#' }
#' Calibration metrics are calculated using each simulation or boot sample. For both options percentile confidence intervals are returned.
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
                          ci = c("sim", "boot", "pw", "none"),
                          n=1000, logitp = T, neval=100, ...){

  # TODO
  # - add ci option "pw" for pointwise (plot only) (DONE)
  # - if smooth = 'none' print metrics for logistic calibration? or remove this option...
  # - vignettes in external and internal validation using pmcalibration
  # - save the boot/sim samples and have summary calculate 95% CIs? (DONE)

  call <- match.call()
  dots <- list(...)

  chk::vld_compatible_lengths(y, p)

  # method <- match.arg(method)
  smooth <- match.arg(smooth)
  ci <- match.arg(ci)

  if (smooth == "none" & isFALSE(logitp)){
    warning("for smooth = 'none' (logistic calibration) logitp is set to TRUE")
    logitp <- T
  }

  if (ci == "boot"){
    if ("cores" %in% names(dots)){
      cores <- dots[["cores"]]
    } else{
      cores <- 1
    }
  }

  # check_method_smooth(method, smooth)

  if (!is.null(neval) && neval > 0){
    pplot <- seq(min(p), max(p), length.out = neval)
  } else{
    pplot <- NULL
  }

  pw <- ci == "pw"
  if (is.null(pplot) & pw){
    warning("ci = 'pw' but neval = 0 or is.null. Will not calculate pointwise standard errors")
    pw <- F
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
    cal <- glm_cal(y = y, p = p, x = x, xp = xp, smooth = smooth,
                   save_data = T, save_mod = T, pw = pw, ...)
  } else if (smooth == "lowess"){
    if (pw) warning("ci = 'pw' is ignored for smooth = 'lowess'")
    cal <- lowess_cal(y = y, p = p, x = x, xp = xp, save_data = T)
  } else if (smooth == "loess"){
    cal <- loess_cal(y = y, p = p, x = x, xp = xp, save_data = T, save_mod = T, pw = pw)
  } else if (smooth == "gam"){
    cal <- gam_cal(y = y, p = p, x = x, xp = xp, save_data = T, save_mod = T, pw = pw, ...)
  }

  if (ci == "boot"){
    b.cal <- run_boots(cal, R = n, cores = cores)
    #conf.int <- get_ci(b.cal, conf_level = conf_level)
  } else if (ci == "sim"){
    b.cal <- simb(cal, R = n)
    #conf.int <- get_ci(b.cal, conf_level = conf_level)
  } else{
    #conf.int <- NULL
    b.cal <- NULL
  }

  if (!is.null(b.cal)){
    #extract samples
    metrics.samples <- do.call(rbind, lapply(b.cal, function(samp) samp$metrics))
    if (any(is.na(metrics.samples))){
      warning("Metrics samples contain NAs. This is probably due to a loess not extrapolating for to-be-plotted values...")
    }

    if (!is.null(b.cal[[1]]$p_c_plot)){
      plot.samples <- do.call(rbind, lapply(b.cal, function(samp) samp$p_c_plot))
      if (any(is.na(plot.samples))){
        warning("Plot samples contain NAs. This is probably due to a loess not extrapolating for to-be-plotted values...")
      }
    } else{
      plot.samples <- NULL
    }
  } else{
    metrics.samples <- NULL
    plot.samples <- NULL
  }

  out <- list(
    call = call,
    metrics = cal$metrics,
    #conf.int = conf.int$metrics,
    metrics.samples = metrics.samples,
    plot = list(
      p = pplot,
      p_c_plot = cal$p_c_plot,
      p_c_plot_se = cal$p_c_plot_se,
      plot.samples = plot.samples
      #conf.int = conf.int$p_c_plot
    ),
    smooth = smooth,
    ci = ci, #conf_level = conf_level,
    n = n,
    smooth_args = cal$smooth_args,
    logitp = logitp
  )

  class(out) <- "pmcalibration"

  return(out)

}
