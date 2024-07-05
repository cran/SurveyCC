#' Plot cc object
#'
#' @param x the survey cc object, produced by [surveycc()]
#' @param dim1,dim2 determines which canonical variates serve as the horizontal
#'  and vertical axes in the plot. Default is dim1 = 1, dim2 = 2. Must not
#'  exceed total number of canonical variates in survey cc object.
#' @param ... Not used.
#'
#' @return the plots
#'
#' @method plot surveycc
#' @export
#' @examples
#' # PATH example
#' design_object <-
#'  survey::svrepdesign(
#'  id = ~PERSONID,
#'  weights = ~R01_A_PWGT,
#'  repweights = "R01_A_PWGT[1-9]+",
#'  type = "Fay",
#'  rho = 0.3,
#'  data=reducedPATHdata,
#'  mse = TRUE
#'  )
#' var.x <- c("R01_AC1022", "R01_AE1022", "R01_AG1022CG")
#' var.y <- c("R01_AX0075", "R01_AX0076")
#' howmany <- 2
#' out <- surveycc(design_object, var.x, var.y, howmany = howmany,
#'   selection = "ROWS")
#' plot(out, dim1 = 1, dim2 = 2)
#'
plot.surveycc <- function(x, dim1 = 1, dim2 = 2, ...) {
  cc_object <- x$cc_object
  #input check
  stopifnot(
    "dim1 & 2 must be positive whole number" = dim1%%1 == 0 && 0 < dim1 && !is.null(dim1),
    "dim1 & 2 must be positive whole number" = dim2%%1 == 0 && 0 < dim2 && !is.null(dim2),
    "dim1 & 2 cannot exceed num canonical variates" = max(dim1, dim2) <= length(x$Stats.cancor),
    "dim1 and dim2 must be different" = dim1 != dim2
  )
  # comment
  s1 <- diag(sqrt(diag(stats::cov(cc_object$X))))
  mycoef1 <- s1 %*% cc_object$coef$X
  stdU <- cc_object$X %*% mycoef1
  s2 <- diag(sqrt(diag(stats::cov(cc_object$Y))))
  mycoef2 <- s2 %*% cc_object$coef$Y
  stdV <- cc_object$Y %*% mycoef2

  graphcoef <- rbind(mycoef1, mycoef2)

  

  
  graphcoef <- cbind(graphcoef, c(rep(0, nrow(mycoef1)), rep(1, nrow(mycoef2))))

  # plot a unit circle
  x <- seq(-1, 1, length = 100)
  y <- sqrt(1 - x^2)
  plot(x, y, type = "l", xlab = paste0("CC", dim1), ylab = paste0("CC", dim2), lty = 2, xlim = c(-2, 2), ylim = c(-2, 2))
  graphics::lines(x, -y, lty = 2)

  graphics::points(
    graphcoef[, dim1], graphcoef[, dim2],
    col = ifelse(graphcoef[, ncol(graphcoef)] == 0, "red", "blue"),
    pch = ifelse(graphcoef[, ncol(graphcoef)] == 0, 1, 3)
  )
  graphics::text(
    graphcoef[, dim1], graphcoef[, dim2],
    labels = c(cc_object$names$X, cc_object$names$Y),
    cex = 0.5, offset = 0.5, pos = 3
  )

  stdV_dim1 <- stdV[, dim1]
  stdU_dim1 <- stdU[, dim1]

  plot(stdU_dim1, stdV_dim1, xlab = paste0("U", dim1), ylab = paste0("V", dim1))
  graphics::text(stdU_dim1, stdV_dim1, labels = seq_along(stdU_dim1))
}
