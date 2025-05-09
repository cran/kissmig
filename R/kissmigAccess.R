#' Get accessiblity map from kissmig output
#'
#' \command{kissmigAccess} calculates an accessibility map from a \code{\link{kissmig}} output of first occurrence (\code{type}=`"FOC"`).
#' These maps allows the integration of limited migration in species distribution models and macroecological analyses.
#' @param grd SpatRaster with one layer as first occurrence generated by kissmig.
#' @param rel bool. If `TRUE`, kissmigAccess returns relative values with maximum `1`, otherwise absolute integer values.
#' Defaults to `FALSE`.
#'
#' @details
#' \code{\link{kissmig}} maps of first occurrences show values of the first iteration step a raster cell was colonized. Early
#' colonized cells have low values, late colonized cells high values. These values are the opposite of accessibility,
#' which is high for early colonized, and low for late colonized cells. \command{kissmigAccess} simply calculates for each
#' cell the accessibility as the difference between the cell value and \command{max(grd)+1}. Cells which have never been
#' colonized remain unchanged (value `0`).
#' @seealso \code{\link{kissmig}}
#' @import terra
#' @export kissmigAccess

kissmigAccess <- function(grd, rel = FALSE) {
  if (!is(grd, "SpatRaster")) stop("parameter 'grd' must be a SpatRaster of the terra package")
  if (dim(grd)[3]!=1) stop("parameter 'grd' must be a SpatRaster with a single layer")
  v <- values(grd)
  v[is.na(v)] <- 0
  values(grd) <- v

  values(grd)[values(grd) > 0] <- max(values(grd)) - values(grd)[values(grd) > 0] + 1

  if (rel == TRUE) {
    grd <- grd / max(values(grd))
  }

  return(grd)
}
