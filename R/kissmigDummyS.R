#' Get a simple suitability map
#'
#' \command{kissmigDummyS} is a support function to generate simple suitability maps based on mean annual air temperature for example code.
#' @usage
#' kissmigDummyS(mean, sd)
#' @param mean temperature mean (degree celsius) of the suitability distribution
#' @param sd temperature standard deviation (degree celsius) of the  suitability distribution
#'
#' @details
#' \command{kissmigDummyS} returns a suitability map for a given extent based on mean annual temperature. It uses
#' data of WorldClim and calculates suitability as a normal distribution defined by \command{mean} and \command{sd}
#' of mean annual temperature. The density function is linearly rescaled to a maximum of 1, the maximum suitability
#' used in \command{kissmig}.
#' @seealso \code{\link{kissmig}}
#' @importFrom stats dnorm
#' @import raster
#' @export kissmigDummyS
#' @references \url{https://www.worldclim.org/}

kissmigDummyS <- function(mean, sd) {
  ans <- wcl_bio1_europe
  if (((mean)<minValue(ans)) | ((mean)>maxValue(ans))) {
    warning("mean outside temperature range")
  }
  values(ans) <- dnorm(values(ans), mean, sd)/dnorm(mean, mean, sd)
  return(ans)
}

utils::globalVariables("wcl_bio1_europe")
