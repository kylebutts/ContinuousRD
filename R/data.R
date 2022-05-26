#' Simulated data for Continuous RD method
#'
#' The simulated data has a true treatment effect of 1 for units below the 30th
#'   percentile. The estimator should be able to detect a treatment effect for
#'   units below the 30th percentile and estimates should be about 1. Note that
#'   since the effect on T of crossing the interval gets smaller as you approach
#'   the 30th percentile, the treatment
#'
#' @format A data frame with 10000 rows and 6 variables:
#' \describe{
#'   \item{id}{Unit identifier}
#'   \item{u}{"unobserved" quantile of T for that unit}
#'   \item{R}{Running variable. Cutoff is at 0}
#'   \item{Z}{Indicator variable for R > 0}
#'   \item{T}{Continuous treatment variable}
#'   \item{Y}{Outcome variable of interest}
#' }
#' @source \url{http://www.diamondse.info/}
"sim_data"
