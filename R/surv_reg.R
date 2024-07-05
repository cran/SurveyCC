#' create survey objects
#'
#' @param n_cc number of canonical correlations
#' @param OGStataUVW an internal object containing the canonical variates
#' @param design_object design object
#'
#' @return a list with the following elements: svyreg1, svyreg2.
#' @noRd
create_survey_objects <- function(n_cc, OGStataUVW, design_object) {

  names_OGStataUVW <- names(OGStataUVW)
  secondindex <- 1:n_cc + n_cc


  # NOTE: we need to modify the design_object to contain the cancor variables
  design_object$variables <- cbind(
    design_object$variables, OGStataUVW[1:n_cc]
  )
  design_object$variables <- cbind(
    design_object$variables, OGStataUVW[secondindex]
  )
  # first regression
  texts1 <- paste(
    names_OGStataUVW[1:n_cc],
    "~", names_OGStataUVW[secondindex]
  )
  formulas1 <- sapply(texts1, stats::as.formula)
  svyreg1 <- lapply(
    formulas1,
    function(formula) survey::svyglm(formula, design = design_object)
  )

  # second regression
  texts2 <- paste(
    names_OGStataUVW[secondindex],
    "~", names_OGStataUVW[1:n_cc]
  )
  formulas2 <- sapply(texts2, stats::as.formula)
  svyreg2 <- lapply(
    formulas2,
    function(formula) survey::svyglm(formula, design = design_object)
  )

  return(list(svyreg1, svyreg2))
}

#' Use the new method to calculate p-values
#'
#' @param n_cc number of canonical correlations
#' @param OGStataUVW data frame with the original data
#' @param design_object the complex survey design object
#'
#' @return a matrix with the following columns: stats, df1, df2, fstat, pval
#' @noRd
weighted_reg <- function(n_cc, OGStataUVW, design_object) {
  
  # creating survey objects
  svyregs <- create_survey_objects(
    n_cc, OGStataUVW,
    design_object
  )
  p1 <- sapply(svyregs[[1]], function(x) stats::coef(summary(x))[2, 4])
  p2 <- sapply(svyregs[[2]], function(x) stats::coef(summary(x))[2, 4])
  svyreg <- lapply(
    1:n_cc,
    function(i) {
      ifelse(p1[i] >= p2[i], svyregs[[1]][i], svyregs[[2]][i])
    }
  )

  # calculating the statistics and p-values
  stats <- sapply(
    1:n_cc,
    function(i) stats::coef(summary(svyreg[[i]][[1]]))[2, 1]
  )
  df1 <- sapply(
    1:n_cc,
    function(i) {
      svyreg[[i]][[1]]$df.null
    }
  )
  df2 <- sapply(
    1:n_cc,
    function(i) svyreg[[i]][[1]]$df.residual
  )
  fstat <- sapply(
    1:n_cc,
    function(i) {
      aux <- summary(svyreg[[i]][[1]])$coefficients[, 3]
      aux[2]
    }
  )
  pval <- sapply(
    1:n_cc,
    function(i) summary(svyreg[[i]][[1]])$coef[2, 4]
  )

  return(cbind(stats, df1, df2, fstat, pval))
}
