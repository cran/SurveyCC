#' Plot canon
#'
#' @param cc_object the cancor object
#' @param dim1 the first dimension
#' @param dim2 the second dimension
#'
#' @return does not return object at present, only plots
#'
#' @noRd
plot_canon <- function(cc_object, dim1, dim2) {
  s1 <- diag(sqrt(diag(stats::cov(cc_object$X))))
  mycoef1 <- s1 %*% cc_object$coef$X
  stdU <- cc_object$X %*% mycoef1
  s2 <- diag(sqrt(diag(stats::cov(cc_object$Y))))
  mycoef2 <- s2 %*% cc_object$coef$Y
  stdV <- cc_object$Y %*% mycoef2

  graphcoef <- rbind(mycoef1, mycoef2)

  #print(graphcoef) #AARON unsure if this is still desired

  #graphcoef <- cbind(graphcoef, c(rep(0, ncol(mycoef1)), rep(1, ncol(mycoef2))))
  graphcoef <- cbind(graphcoef, c(rep(0, nrow(mycoef1)), rep(1, nrow(mycoef2))))

  # grid for two plots
  #graphics::par(mfrow = c(1, 2)) #2023-09-08 cancelling par manipulation as per CRAN email

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
    labels = c(cc_object$names$X, cc_object$names$Y)
  )

  stdV_dim1 <- stdV[, dim1]
  stdU_dim1 <- stdU[, dim1]

  plot(stdU_dim1, stdV_dim1, xlab = paste0("U", dim1), ylab = paste0("V", dim1))
  graphics::text(stdU_dim1, stdV_dim1, labels = 1:length(stdU_dim1))

  #graphics::par(mfrow = c(1, 1))
}
