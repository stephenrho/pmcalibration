#' Simulation based inference with a calibration curve object
#'
#' @param x an object created using one of the \code{cal} functions
#' @param ... other object specific arguments
#'
#' @return simulated calibration metrics and values for plotting
#' @keywords internal
#' @export
simb <- function(x, ...){
  UseMethod("simb", x)
}

#' @rdname simb
#' @export
simb.glm_cal <- function(x, R = 1000){

  mod <- x$model; y = x$y; p = x$p; X = x$X; Xp = x$Xp

  b <- coef(mod); Vb <- vcov(mod)

  B <- MASS::mvrnorm(n = R, mu = b, Sigma = Vb)

  out <- lapply(1:R, function(i){
    p_c <- as.numeric(invlogit( as.matrix(cbind(1, X)) %*% B[i, ]))

    if (!is.null(Xp)){
      p_c_plot <- as.numeric(invlogit( as.matrix(cbind(1, Xp)) %*% B[i, ]))
    } else{
      p_c_plot <- NULL
    }

    s <- list(
      p_c = p_c,
      metrics = cal_metrics(p, p_c),
      p_c_plot = p_c_plot
    )

    return(s)
  }
  )

  return(out)
}

#' @rdname simb
#' @export
simb.lowess_cal <- function(x, ...){
  stop("ci = 'sim' is not available for lowess or loess calibration curves (try ci = 'boot')")
}

#' @rdname simb
#' @export
simb.loess_cal <- function(x, ...){
  stop("ci = 'sim' is not available for lowess or loess calibration curves (try ci = 'boot')")
}

#' @rdname simb
#' @export
simb.gam_cal <- function(x, ...){
  x
}
