#' Canonical correlation analysis for complex survey data
#'
#' @description This command extends the functionality of [candisc::cancor] by
#'  calculating the test statistics, degrees of freedom and p-values necessary
#'  to estimate and interpret the statistical significance of the secondary canonical
#'  corr according to the methods Wilks' lambda, Pillai's trace, and Hotelling-Lawley
#'  trace (Caliński et al., 2006) and Roy's largest root (Johnstone, 2009). The units
#'  and variables graphs (Gittins, 1986) can also be drawn by `surveycc` further
#'  complementing the information listed by the existing `cancor`.
#'
#'  Moreover, `csdcanon` implements an algorithm (Cruz-Cano, Cohen, and Mead-Morse, 2024) that
#'  allows the inclusion of complex survey design elements, e.g. strata, cluster and
#'  replicate weights, in the estimation of the statistical significance of the canonical
#'  correlations. The core idea of the algorithm is to reduce the problem of finding
#'  the correlations among the canonical variates and their corresponding statistical
#'  significance to calculating an equivalent sequence of univariate linear regression.
#'  This switch allows the user to take advantage of the existing theoretical and
#'  computational resources that integrate the complex survey design elements into
#'  these regression models (Valliant and Dever, 2018). Hence, this algorithm can include
#'  the same complex design elements as in `survey`.
#'
#' @param design_object a survey design object generated from package `survey`,
#'  eg [survey::svydesign]
#' @param var.x the first set of variables; a vector of names
#' @param var.y the second set of variables; a vector of names
#' @param howmany positive integer; allows the user to choose the number of canonical correlations
#'  for which the statistical significance statistics are displayed. Default is to choose
#'  the minimum of `length(var.x)` and `length(var.y)`. Cannot exceed this value.
#' @param dim1,dim2 determines which canonical variates serve as the horizontal
#'  and vertical axes in the optional plot. NOTE: if dim1 and dim2 not provided,
#'  no graph will be displayed.
#' @param selection allows the user to choose whether the type of sample size is equal to
#'  the number of rows in the data set or the sum of the survey weights.
#'
#' @return A list, containing the canonical correlation object, as well as tables
#'  of the various tests of significance. This includes the test statistics,
#'  degrees of freedom, and p-values for:
#'  * Wilk's lambda
#'  * Pillai's trace
#'  * Hotelling-Lawley
#'  * Roy's greatest root
#'  * the Cruz-Cano algorithm using the survey design object
#'
#'  NOTE: For more information on the statistics presented, i.e. test statistic, df1,
#'  df2, Chi-Sq/F and p-val, please see the documentation in [candisc::cancor]
#'  for Wilk's Lambda, Pillai's Trace and Hotelling-Lawley Trace
#'  (although the present package uses a Chi-squared approximation to the F-distribution),
#'  and see the documentation in [survey::svyglm] for the Weighted/Complex Survey Design
#'  regression.
#'
#' @references  * Cruz-Cano, Cohen, and Mead-Morse. Canonical Correlation Analysis of
#'  Survey data: The SurveyCC R package. The R Journal under review; 2024.
#'  * Gentzke AS, Wang TW, Cornelius M, Park-Lee E, Ren C, Sawdey MD, Cullen KA,
#'  Loretan C, Jamal A, Homa DM. Tobacco Product Use and Associated Factors among
#'  Middle and High School Students - National Youth Tobacco Survey, United States, 2021.
#'  rveill Summ. 2022;71(5):1-29. doi: 10.15585/mmwr.ss7105a1. PubMed PMID: 35271557.
#'  * Gittins R.  Canonical Analysis: A Review with Applications in Ecology:
#'  Springer Berlin Heidelberg; 1986.
#'  * Caliński T., Krzyśko M. and WOłyński W. (2006) A Comparison of Some Tests
#'  for Determining the Number of Nonzero Canonical Correlations, Communications
#'  in Statistics -Simulation and Computation, 35:3, 727-749, DOI: 10.1080/03610 6290.
#'  * Hyland A, Ambrose BK, Conway KP, et al. Design and methods of the Population
#'  Assessment of Tobacco and Health (PATH) StudyTobacco Control 2017;26:371-378.
#'  * Johnstone IM. Approximate Null Distribution of the largest root in a
#'  Multivariate Analysis.  Ann Appl Stat. 2009;3(4):1616-1633. doi: 10.1214/08-AOAS220.
#'  PMID: 20526465; PMCID: PMC2880335.
#'  *  Valliant R. and Dever JA.  Survey Weights: A Step-by-Step Guide to Calculation:
#'  Stata Press; 2018. ISBN-13: 978-1-59718-260-7.
#'
#' @export
#'
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
#' dim1 <- 1
#' dim2 <- 2
#' surveycc(design_object, var.x, var.y, howmany = howmany,
#'   dim1 = dim1, dim2 = dim2, selection = "x")
#'
#' # NYTS example
#' design_object <-
#'   survey::svydesign(
#'   ids = ~psu2,
#'   nest = TRUE,
#'   strata = ~v_stratum2,
#'   weights = ~finwgt,
#'   data = reducedNYTS2021data
#' )
#' var.x <- c("qn9", "qn38", "qn40", "qn53", "qn54", "qn64", "qn69", "qn74",
#'            "qn76", "qn78", "qn80", "qn82", "qn85", "qn88", "qn89")
#' var.y <- c("qn128", "qn129", "qn130", "qn131", "qn132", "qn134")
#' howmany <- 3
#' surveycc(design_object = design_object, var.x = var.x,
#'   var.y = var.y, howmany = howmany, selection = "x")
#'
surveycc <- function(
    design_object, var.x, var.y, howmany = NA, dim1 = NA, dim2 = NA,
    selection = "FREQ") {
  # modify the design object to have complete cases on var.x and var.y
  design_object <- subset(
    design_object,
    stats::complete.cases(
      design_object$variables[, c(var.x, var.y)]
    )
  )
  # This is an attempt to extract weights for the cancor object
  if (!is.null(design_object$pweights)) {
    myweights <- design_object$pweights
  } else myweights <- weights(design_object)

  # create a version of the svy-obj w/ just weights, for weighted regression
  design_objectW <- survey::svydesign(
    ids = ~1, data = design_object$variables, weights = myweights
  )
  # create the cancor object
  cc_object <- candisc::cancor(
    design_object$variables[, var.x],
    design_object$variables[, var.y],
    weights = myweights
  )

  # draw plot
  # NOTE AARON 2023-08-09 only plot if dim1 and dim2 provided
  if (!is.na(dim1) & !is.na(dim2)) {
    plot_canon(cc_object, dim1, dim2)
  }

  # preparing the inputs for the calcpval function
  varcount1 <- length(cc_object$names$X)
  varcount2 <- length(cc_object$names$Y)
  if (varcount1 >= varcount2) {
    OgX <- cc_object$X
    OgY <- cc_object$Y
    myrawcoef_var1 <- cc_object$coef$X
    myrawcoef_var2 <- cc_object$coef$Y
  } else {
    OgX <- cc_object$Y
    OgY <- cc_object$X
    myrawcoef_var2 <- cc_object$coef$X
    myrawcoef_var1 <- cc_object$coef$Y
  }
  # if (command != "classic") { # no longer necessary, DELETE
  #   diag_W <- diag(cc_object$weights)
  # } else {
  #   diag_W <- diag(nrow = nrow(OgX))
  # }
  diag_W <- diag(cc_object$weights)
  n_cc <- cc_object$ndim
  weights <- cc_object$weights

  # calculating and formating the p-values
  output <- calcpval(
    cc_object, selection, OgX, OgY, diag_W,
    myrawcoef_var1, myrawcoef_var2, weights
  )
  Results <- output$Results
  OGStataUVW <- output$OGStataUVW


  # applying the new method
  w_reg_results <- weighted_reg(n_cc, OGStataUVW, design_objectW)
  svy_reg_results <- weighted_reg(n_cc, OGStataUVW, design_object)

  # adding it to results
  for (i in 1:n_cc) {
    Results[[i]] <- rbind(Results[[i]], w_reg_results[i, ])
    Results[[i]] <- rbind(Results[[i]], svy_reg_results[i, ])
    rownames(Results[[i]]) <- c(
      rownames(Results[[i]])[1:4],
      "Weighted Survey CC",
      "Complex Survey CC"
    )
  }

  if (is.na(howmany)) {
    howmany <- n_cc
  }

  names(Results) <- paste(
    "Stats.cancor.", 1:n_cc,
    sep = ""
  )

  # NOTE 2023-08-20 round digits for tables
  Results <- lapply(Results, round, digits = 5)

  # printing the results
  # for (i in 1:howmany) {
  #   aux <- Results[[i]][[1]]
  #   mag <- round(aux[1], 0.0001)
  #   ResultsAux <- Results[[i]]
  #   print(ResultsAux,
  #         title = paste0("Statistics for Canonical Correlation: ", i),
  #         twidth = 30, format = "%10.4f",
  #         rowtitle = paste0("Canonical Correlation=", mag),
  #         border = c("top", "bottom")
  #   )
  # }

  out.results <- list(
    Stats.cancor = Results[1:howmany],
    cc_object = cc_object
  )

  return(out.results)
}
