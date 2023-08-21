#' Create a calibration curve
#'
#' @description
#' Assess calibration of clinical prediction models (agreement between predicted and observed probabilities) via different smooths.
#'
#' @param y a binary or time-to-event outcome. Latter must be an object created via \code{survival::Surv}. Only right censored outcomes.
#' @param p predicted probabilities from a clinical prediction model. For a time-to-event object \code{time} must be specified and \code{p} are predicted probabilities of the outcome happening by \code{time} units of time follow-up.
#' @param smooth what smooth to use. Available options:
#' \itemize{
#' \item{'rcs' = restricted cubic spline using \code{rms::rcs}. Optional arguments for this smooth are \code{nk} (number of knots; defaults to 5) and \code{knots} (knot positions; set by \code{Hmisc::rcs.eval} if not specified) }
#' \item{'ns' = natural spline using \code{splines::ns}. Optional arguments are \code{df} (default = 6), \code{knots}, \code{Boundary.knots} (see \code{?splines::ns}) }
#' \item{'bs' = B-spline using \code{splines::bs}. Optional arguments are \code{df} (default = 6), \code{knots}, \code{Boundary.knots} (see \code{?splines::bs}) }
#' \item{'gam' = generalized additive model via \code{mgcv::gam} and \code{mgcv::s}. Optional arguments are \code{bs}, \code{k}, \code{fx}, \code{method} (see \code{?mgcv::gam} and  \code{?mgcv::s}) }
#' \item{'lowess' = uses \code{lowess(x, y, iter = 0)} based on \code{rms::calibrate}. Only for binary outcomes.}
#' \item{'loess' = uses \code{loess} with all defaults. Only for binary outcomes. }
#' \item{'none' = logistic regression with single predictor variable (performs logistic calibration when \code{logitp = T}). See \code{\link{logistic_cal}} }
#' }
#' 'rcs', 'ns', 'bs', and 'none' are fit via \code{glm} or \code{survival::coxph} and 'gam' is fit via \code{mgcv::gam} with \code{family = Binomial(link="logit")} for a binary outcome or \code{mgcv::cox.ph} when \code{y} is time-to-event.
#' @param time what follow up time do the predicted probabilities correspond to? Only used if \code{y} is a \code{Surv} object
#' @param ci what kind of confidence intervals to compute?
#' \itemize{
#' \item{'sim' = simulation based inference; \code{n} samples are taken from a multivariate normal distribution with mean vector = coef(mod) and variance covariance = vcov(model).}
#' \item{'boot' = bootstrap resampling with \code{n} replicates. \code{y} and \code{p} are sampled with replacement and calibration curve is reestimated. If \code{knots} are specified the same knots are used for each resample (otherwise they are calculated using resampled x)}
#' \item{'pw' = pointwise confidence intervals calculated via the standard errors produced by relevant \code{predict} methods. Only for plotting curves - CIs are not produced for metrics (not available for smooth = 'lowess')}
#' }
#' Calibration metrics are calculated using each simulation or boot sample. For both options percentile confidence intervals are returned.
#' @param n number of simulations or bootstrap resamples
#' @param transf transformation to be applied to \code{p} prior to fitting calibration curve. Valid options are 'logit', 'log-log', 'none', or a function (must retain order of \code{p}). If unspecified defaults to 'logit' for binary outcomes and 'log-log' for time-to-event outcomes.
#' @param neval number of points (equally spaced between \code{min(p)} and \code{max(p)}) to evaluate for plotting (0 or NULL = no plotting). Can be a vector of probabilities.
#' @param ... additional arguments for particular smooths. For ci = 'boot' the user is able to run samples in parallel (using the \code{parallel} package) by specifying a \code{cores} argument
#'
#' @references Austin P. C., Steyerberg E. W. (2019) The Integrated Calibration Index (ICI) and related metrics for quantifying the calibration of logistic regression models. \emph{Statistics in Medicine}. 38, pp. 1–15. https://doi.org/10.1002/sim.8281
#' @references Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina M., Steyerberg E.W. (2016). A calibration hierarchy for risk models was defined: from utopia to empirical data. \emph{Journal of Clinical Epidemiology}, 74, pp. 167-176
#' @references Austin, P. C., Harrell Jr, F. E., & van Klaveren, D. (2020). Graphical calibration curves and the integrated calibration index (ICI) for survival models. \emph{Statistics in Medicine}, 39(21), 2714-2742.
#'
#' @returns a \code{pmcalibration} object containing calibration metrics and values for plotting
#' @export
pmcalibration <- function(y, p,
                          smooth=c("none", "ns", "bs", "rcs", "gam", "lowess", "loess"),
                          time = NULL,
                          ci = c("sim", "boot", "pw", "none"),
                          n=1000,
                          transf = NULL,
                          neval=100, ...){

  # TODO
  # - survival methods (DONE glm; todo - gam, hare, flexsurv?)
  # - update summary/print methods for tte outcomes...

  call <- match.call()
  dots <- list(...)

  chk::vld_compatible_lengths(y, p)

  # method <- match.arg(method)
  smooth <- match.arg(smooth)
  ci <- match.arg(ci)

  # if (smooth == "none" & isFALSE(logitp)){
  #   warning("for smooth = 'none' (logistic calibration) logitp is set to TRUE")
  #   logitp <- T
  # }

  surv <- is(y, "Surv")

  if (surv & ci == "sim") stop("Simulation based inference not available for time-to-event outcomes")
  if (surv & smooth %in% c("lowess", "loess")) stop("smooth = 'lowess' or 'loess' not available for time-to-event outcomes")
  if (surv & smooth == "gam") stop("gam currently not implemented for Surv outcomes")

  if (ci == "boot"){
    if ("cores" %in% names(dots)){
      cores <- dots[["cores"]]
    } else{
      cores <- 1
    }
  }

  # check_method_smooth(method, smooth)
  if (length(neval) > 1){
    if (any(neval < 0) | any(neval > 1)){
      stop("When providing p for plotting directly via neval, all elements should be between 0 and 1")
    }
    pplot <- neval
  } else{
    if (!is.null(neval) && neval > 0){
      pplot <- seq(min(p), max(p), length.out = neval)
    } else{
      pplot <- NULL
    }
  }

  pw <- ci == "pw"
  if (is.null(pplot) & pw){
    warning("ci = 'pw' but neval = 0 or is.null. Will not calculate pointwise standard errors")
    pw <- F
  }

  # if (logitp){
  #   xp <- if (!is.null(pplot)) logit(pplot) else NULL
  #   x <- logit(p)
  # } else{
  #   x <- p; xp <- pplot
  # }
  if (is.null(transf)){
    if (surv){
      transf <- "log-log"
    } else{
      transf <- "logit"
    }
  }

  if (is.function(transf)){
    if (length(formals) > 1) stop("transf should just take one argument")
    tfun <- transf
  } else{
    if (!transf %in% c("none", "logit", "log-log")) {
      stop("invalid transf. Should be 'none', 'logit', 'log-log' or a function")
    }
    if (transf == "none"){
      tfun <- function(x) x
    } else if (transf == "logit"){
      tfun <- logit
    } else if (transf == "log-log"){
      tfun <- function(x) log(-log(1 - x))
    }
  }
  xp <- if (!is.null(pplot)) tfun(pplot) else NULL
  x <- tfun(p)

  # fit cc
  if (smooth %in% c("none", "ns", "bs", "rcs")){
    cal <- glm_cal(y = y, p = p, x = x, xp = xp, smooth = smooth, time=time,
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
    transf = transf,
    time = time,
    outcome = cal$outcome
  )

  class(out) <- "pmcalibration"

  return(out)

}
