#' Define a geographic origin
#'
#' \command{kissmigOrigin} is a support function to define a rectangular origin for a \code{\link{kissmig}} call.
#' @param grd SpatRaster with one layer as reference.
#' @param x integer as lower left x-coordinate of geographic origin.
#' @param y integer as lower left y-coordinate of geographic origin.
#' @param size number as size of the quadratic origin.
#'
#' @details
#' \code{\link{kissmigOrigin}} returns a SpatRaster with one layer characterized by the reference
#' \code{grd}. Cell values are set to `0`, except for cells of the origin defined by
#' \code{\link[terra]{ext}(x, x+size, y, y+size)} which are set to `1`.
#' @seealso \code{\link{kissmig}}
#' @import terra
#' @export kissmigOrigin

kissmigOrigin <- function(grd, x, y, size) {
  if (!is(grd, "SpatRaster")) stop("the reference 'grd' must be a SpatRaster of the terra package")
  if (dim(grd)[3]!=1) stop("the reference 'grd' must be a SpatRaster with a single layer")
  ans <- grd
  values(ans) <- 0
  values(ans)[cells(ans, ext(x, x + size, y, y + size))] <- 1.0
  return(ans)
}
