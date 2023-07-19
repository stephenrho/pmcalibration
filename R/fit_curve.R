
fit_curve <- function(y, p, logit, method, usep = F, df = NULL, knots = NULL, k = NULL){
  if (missing(p) & missing(logit)){
    stop("One of p or logit must be specified.")
  }

  if (missing(logit) & !usep){
    logit <- log(p/(1 - p))
  }

  if (missing(p) & usep){
    p <- 1/(1 + exp(logit)) # check this!!
  }

  if (length(method) > 1){
    method <- method[1]
  }
  stopifnot(method %in% c("natural", "cubic_reg", "cubic_pen", "tp", "lowess", "loess"))

  if (usep){
    x <- p
  } else{
    x <- logit
  }

  if (method == "natural"){
    if (is.null(df) && is.null(knots)){
      df <- 5
    }
    X <- splines::ns(x, df = df, knots = knots)

    kn <- c(attributes(X)$Boundary.knots[1], attributes(X)$knots, attributes(X)$Boundary.knots[2])
    cat("Fitting a calibration curve on",
        ifelse(usep, "p", "logit(p)"),
        "via a natural cubic spline with df =",
        attributes(X)$dim[2],
        "and knots at",
        kn)

    cc <- glm(y ~ X, family = "binomial")
  }
  if (method %in% c("cubic_reg", "cubic_pen", "tp")){
    if (is.null(k)){
      if (is.null(df)){
        k <- 10
      } else{
        k <- df
      }
    }

    if (method %in% c("cubic_reg", "cubic_pen")){
      bs <- "cr"
    } else if (method == "tp"){
      bs <- "tp"
    }

    cat("Fitting a calibration curve on",
        ifelse(usep, "p", "logit(p)"),
        "via a",
        if (method == "cubic_reg") "cubic regression" else if (method == "cubic_pen") "penalized cubic" else if (method == "tp") "thin plate regression",
        "spline with k =",
        k)

    cc <- mgcv::gam(y ~ s(x, k = k, fx = method == "cubic_reg", bs = bs), family = "binomial")
  }
  return(cc)
}
