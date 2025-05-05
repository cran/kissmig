#' Get a simple suitability map
#'
#' \command{kissmigDummyS} is a support function to generate simple climate suitability maps
#' for example code.
#' @param mean number of the mean temperature (degree celsius) of the suitability distribution.
#' @param sd number of the temperature standard deviation (degree celsius) of the  suitability distribution.
#'
#' @details
#' \command{kissmigDummyS} returns a SpatRaster of a simple climate suitability map for Europe based on mean annual temperature. It uses
#' data of WorldClim and calculates suitability as a normal distribution defined by \code{mean} and \code{sd}
#' of mean annual temperature. The density function is linearly rescaled to a maximum of `1`, the maximum suitability
#' used in \code{\link{kissmig}}.
#' @seealso \code{\link{kissmig}}
#' @importFrom stats dnorm
#' @import terra
#' @export kissmigDummyS
#' @references \url{https://www.worldclim.org/}

kissmigDummyS <- function(mean, sd) {
  ans <- terra::rast(kissmig::wcl_bio1_europe)

  if (((mean) < minmax(ans)[1]) | ((mean) > minmax(ans)[2])) {
    warning("mean outside temperature range")
  }

  values(ans) <- dnorm(values(ans), mean, sd) / dnorm(mean, mean, sd)

  return(ans)
}
