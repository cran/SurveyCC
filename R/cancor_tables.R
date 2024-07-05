#' computes a normalized version of the inner product defined
#' by a symmetric matrix
#'
#' @param x vector
#' @param y vector
#' @param A symmetric matrix
#'
#' @return the normalized inner product
#' @noRd
norm_inner_prod <- function(x, y, A) {
  crossprod(x, A %*% y) / sum(diag(A))
}

#' Computes the weighted rhos
#'
#' @param U matrix of the first set of canonical vectors
#' @param V matrix of the second set of canonical vectors
#' @param diag_W diagonal matrix of the weights
#'
#' @return a vector of the weighted rhos
#' @noRd
calculate_weighted_rho <- function(U, V, diag_W) {
  # for each column of U and V substract its mean
  cen_U <- sweep(U, 2, colMeans(U), "-")
  cen_V <- sweep(V, 2, colMeans(V), "-")

  # compute the weighted variances
  weighted_varU <- diag(norm_inner_prod(cen_U, cen_U, diag_W))
  weighted_varV <- diag(norm_inner_prod(cen_V, cen_V, diag_W))
  denom <- sqrt(weighted_varU * weighted_varV)
  weighted_rho <- diag(norm_inner_prod(cen_U, cen_V, diag_W)) / denom

  return(weighted_rho)
}


#' Computes the Wilks' Lambda test statistic and p-value for each
#' canonical correlation.
#'
#' @param weighted_rho A vector of weighted canonical correlations.
#' @param p The number of variables in the first set.
#' @param q The number of variables in the second set.
#' @param Lawley Lawley's term
#' @param samplesize The sample size.
#'
#' @return A matrix with the following columns: Lambda, df, NA, ChiSq_Wilks_Lambda, p_val_Wilks_Lambda.
#' @noRd
table_wilks <- function(weighted_rho, p, q, samplesize, Lawley) {
  s <- length(weighted_rho)
  Lambda <- rev(cumprod(rev(1 - (weighted_rho^2))))
  df <- (p + 1 - 1:s) * (q + 1 - 1:s)
  ChiSq_Wilks_Lambda <- -((samplesize - 1) - 0.5 * (p + q + 1) + Lawley - 1:s) * log(Lambda)
  p_val_Wilks_Lambda <- stats::pchisq(ChiSq_Wilks_Lambda, df, lower.tail = FALSE)
  return(cbind(Lambda, df, NA, ChiSq_Wilks_Lambda, p_val_Wilks_Lambda))
}

#' compute Roys test statistic and p-value
#'
#' @param weighted_rho vector of weighted canonical correlations
#' @param samplesize effective sample size
#' @param p number of variables in the first set
#' @param q number of variables in the second set
#'
#' @return a matrix with the following columns: largest_root, v1, v2,
#' Roys_Greatest_Root_stat, Roys_Greatest_Root_p_value
#' @noRd
table_roys <- function(weighted_rho, samplesize, q, p) { #NOTE going to try swapping p and q
  # parameters
  s <- length(weighted_rho) - 1
  par1 <- rep(p, s + 1)
  par2 <- samplesize + 0:s - q - 1
  par3 <- q - 0:s
  fraks <- pmin(par1, par3)
  frakm <- (abs(par1 - par3) - 1) / 2
  frakn <- (par2 - par1 - 1) / 2

  # degrees of freedom
  v1 <- fraks + (2 * frakm) + 1
  v2 <- fraks + (2 * frakn) + 1

  # stats and p-values
  largest_root <- rep(weighted_rho[1]^2, s + 1)
  stat <- (v2 * largest_root) / (v1 * (1.0 - largest_root))
  p_value <- stats::pf(stat, v1, v2, lower.tail = FALSE)

  return(cbind(
    largest_root, v1, v2, stat, p_value
  ))
}

#' compute Pillai's trace statistic
#'
#' @param weigthed_rho vector of weighted canonical correlations
#' @param Lawley vector of Lawley's terms
#' @param p number of variables in X
#' @param samplesize The sample size
#' @param q number of variables in Y
#'
#' @return a matrix with the Pillai's trace statistic and its p-value
#' @noRd
table_pillais <- function(weigthed_rho, Lawley, p, q, samplesize) {
  s <- length(weigthed_rho)
  V <- rev(cumsum(rev(weigthed_rho^2)))
  Pillais_Trace_stat <- (samplesize - 1 - 2 * (1:s) + Lawley) * V
  df1 <- (p + 1 - (1:s)) * (q + 1 - (1:s))
  Pillais_Trace_p_value <- stats::pchisq(Pillais_Trace_stat, df1,
                                         lower.tail = FALSE
  )
  return(cbind(V, df1 - 1, NA, Pillais_Trace_stat, Pillais_Trace_p_value))
}

#' calculate the Hotelling-Lawley Trace statistic
#'
#' @param weighted_rho a vector of weighted rhos
#' @param p number of variables in X
#' @param q number of variables in Y
#' @param samplesize effective sample size
#' @param Lawley a vector of Lawley's terms
#' @noRd
table_hotelling <- function(weighted_rho, p, q, samplesize, Lawley) {
  s <- length(weighted_rho)
  U <- rev(cumsum(rev(weighted_rho^2 / (1 - weighted_rho^2))))
  Hotelling_Lawley_Trace_stat <- (samplesize - p - q - 2 + Lawley) * U
  df1 <- (p + 1 - (1:s)) * (q + 1 - (1:s))
  Hotelling_Lawley_Trace_p_value <- stats::pchisq(Hotelling_Lawley_Trace_stat, df1,
                                                  lower.tail = FALSE
  )
  return(cbind(
    U, df1, NA, Hotelling_Lawley_Trace_stat,
    Hotelling_Lawley_Trace_p_value
  ))
}


#' calculate p-values for the canonical correlations
#'
#' @param cc_object output from the candisc function
#' @param selection "FREQ" or "ROWS"
#' @param OgX original X matrix
#' @param OgY original Y matrix
#' @param diag_W diagonal matrix of weights
#' @param myrawcoef_var1 raw canonical coefficients for the first set
#' @param myrawcoef_var2 raw canonical coefficients for the second set
#' @param weights vector of weights
#'
#' @return a list with the following elements: OGStataU, OGStataV,
#' OGStataUVW, Results.
#' @noRd
calcpval <- function(
    cc_object, selection, OgX, OgY, diag_W,
    myrawcoef_var1, myrawcoef_var2, weights) {
  # reserving space for the results
  output <- list()

  # calculating the canonical vectors
  U <- OgX %*% myrawcoef_var1
  V <- OgY %*% myrawcoef_var2

  # storing the canonical vectors
  output$OGStataU <- U
  output$OGStataV <- V
  output$OGStataUVW <- as.data.frame(cbind(U, V, weights))

  # calculating the weighted rhos
  weighted_rho <- calculate_weighted_rho(U, V, diag_W)

  # calculating the effective sample size
  if (selection == "FREQ") {
    samplesize <- sum(trunc(cc_object$weights))
    samplesize
  } else {
    samplesize <- nrow(cc_object$X)
    samplesize
  }

  # calculating the Lawley's term
  Lawley <- cumsum(1 / (weighted_rho^2))

  # calculating the Wilks Lambda
  wilks <- table_wilks(
    weighted_rho, ncol(OgX), ncol(OgY),
    samplesize, Lawley
  )

  # calculating the Pillai's trace
  pillais <- table_pillais(
    weighted_rho, Lawley, ncol(OgX), ncol(OgY), samplesize
  )

  # calculating the Hotelling-Lawley trace
  hotelling <- table_hotelling(
    weighted_rho, ncol(OgX), ncol(OgY),
    samplesize, Lawley
  )

  # calculating the Roy's greatest root
  roys <- table_roys(
    weighted_rho, samplesize, ncol(OgX), ncol(OgY)
  )

  # putting it all together
  Results <- vector("list", length(weighted_rho))
  for (i in seq_along(weighted_rho)) {
    myResults <- matrix(nrow = 0, ncol = 5)
    #myResults <- rbind(myResults, rep(weighted_rho[i]^2, 5)) #AARON 2023-08-17 comment out
    myResults <- rbind(myResults, wilks[i, ])
    myResults <- rbind(myResults, pillais[i, ])
    myResults <- rbind(myResults, hotelling[i, ])
    myResults <- rbind(myResults, roys[i, ])
    colnames(myResults) <- c(
      "Statistic", "df1", "df2",
      "Chi-Sq/F", "p-val"
    )
    rownames(myResults) <- c(
      "Wilks' Lambda", "Pillai's Trace",
      "Hotelling-Lawley Trace", "Roy's Greatest Root"
    )
    Results[[i]] <- myResults
  }

  # returning the results
  output$Results <- Results
  return(output)
}
