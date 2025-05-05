#' Run a simple dynamic species distribution model
#'
#' \command{kissmig} runs a simple, raster-based, stochastic model to simulate species
#' distribution dynamics driven by environmental suitability and the migration ability of the species.
#' @param O SpatRaster with a single layer of the initial species distribution as the geographic origin.
#' @param S SpatRaster with a single or multiple suitability layers. For the default setting `NULL` the suitability of all cells is set to 1.0 (maximum).
#' @param it integer of the number of iteration steps representing the species' migration ability and applied to each suitability layer of S.
#' @param type string defining the type of resulting map. Default setting is `"FOC"`.
#'    * `"DIS"` final species distribution after the last iteration step
#'    * `"FOC"` number of the iteration step of the first occurrence
#'    * `"LOC"` number of the iteration step of the last occurrence
#'    * `"NOC"` number of the iteration steps with occurrence
#' @param signed bool. If `TRUE`, the sign indicates whether the cells was colonized (positive) or uncolonized (negative) after the last iteration step.
#' Used in combination with `"FOC"`, `"LOC"`, or `"NOC"`. Default setting is `FALSE`.
#' @param pext probability \[0,1\] that a colonized cell becomes uncolonized between iteration steps, i.e., without recolonization the species gets locally extinct.
#' Default setting is `1.0`.
#' @param pcor probability \[0,1\] that corner cells are considered in the 3x3 cell neighborhood.
#' Default setting is `0.2`.
#' @param seed integer used to set the seed of the random number generator.
#' Default setting is `NULL`.
#' @param n_threads integer of the number of threads for parallel computing.
#' Default setting is `1` (i.e., no parallel computing).
#' @param n_random integer defining the amount of random numbers for the simulation.
#' Default setting is `10000`.
#'
#' @details
#' Starting from an initial species distribution \code{O} as the geographic origin, \command{kissmig} simulates species distributions
#' in an environment characterized by a single or multiple (time series) suitability layers \code{S}.
#' The simulation uses a simple 3x3 cell algorithm, applied \code{it} iteration steps. Colonized cells
#' have value `1`, uncolonized cells value `0`, and environmental suitability values vary binary or
#' continuously between `0` (unsuitable) and `1` (suitability maximum). In case \code{S} consists of
#' multiple suitability layers, \code{it} is applied to each layer.
#'
#' Between iteration steps, colonized cells become extinct with the probability \code{pext}, and for a recolonization
#' or new colonization event corner cells within the 3x3 neighborhood are considered with probability \code{pcor}
#' (\code{pcor}=`0.2` produces more realistic circular spread patterns in homogeneous environments - see Nobis and Normand 2014, \doi{10.1111/ecog.00930}).
#'
#' To get reproducible results, the seed of the R random number generator can be set using the \code{seed} parameter.
#' In addition, for reproducibility the value of \code{n_random} needs to stay the same between simulations.
#' Changing one of these parameters will generate stochastically different results - even with binary suitability due to the probability \code{pcor}
#' of the corner events.
#'
#' Signed results with \code{signed}=`TRUE` can be generated to get in addition to the result type
#' `"FOC"`, `"LCO"`, or `"NOC"` with the same call the final distribution (`"DIS"`) with positive values being colonized and negative values being
#' previously colonized but uncolonized after the last iteration step.
#'
#' Runtime optimization by parallel computing can be initialized by increasing \code{n_threads}. Depending on the hardware used, it is advantageous to set
#' \code{n_threads} significantly higher than the number of available CPU cores. To figure out optimal settings compare execution
#' time with \code{n_threads} set to 2^x (e.g. 4, 8, 16, 32, ...). The implementation uses the OpenMP specification for parallel programming,
#' \url{https://www.openmp.org/}.
#'
#' @references \itemize{
#'    \item{Nobis, M.P. and Normand, S. 2014. KISSMig - a simple model for R to account for limited migration in
#'    analyses of species distributions. \cite{Ecography} 37: 1282-1287. \doi{10.1111/ecog.00930}}
#'    \item{KISSMig homepage \url{https://purl.oclc.org/wsl/kissmig}.}
#' }
#' @seealso \code{\link{kissmigAccess}, \link{kissmigOrigin}}
#' @examples
#' \donttest{
#' # generate some SpatRaster data to run kissmig
#'
#' s <- kissmigDummyS(mean = 12, sd = 3) # suitability map
#' o <- kissmigOrigin(s, 8.0, 44.5, 0.5) # geographic origin
#' l <- s >= 0 # land mask used for visualization below
#' plot(s, asp = 1.0, main = "Climate suitability & origin (in red)")
#' plot(o, col = c(rgb(0,0,0,0), "red"), legend = FALSE, add = TRUE) # add origin
#'
#' # run kissmig with different types of output
#'
#' k <- kissmig(o, s, it = 150, type = "DIS")
#' plot(k*l, asp = 1.0, main = "Final distribution (DIS)")
#' plot(o, col = c(rgb(0,0,0,0), "red"), legend = FALSE, add = TRUE) # add origin
#'
#' sb <- s>0.5 # binary suitability for an arbitrary threshold of 0.5
#' k <- kissmig(o, sb, it = 150, type = "DIS")
#' plot(k*l, asp = 1.0, main = "Final distribution (DIS) using binary suitability")
#' plot(o, col = c(rgb(0,0,0,0), "red"), legend = FALSE, add = TRUE) # add origin
#'
#' k <- kissmig(o, s, it = 150, type = "FOC")
#' plot(k*l, asp = 1.0, main = "First iteration step of occurrence (FOC)",
#'      col = c("lightgrey", map.pal("viridis", n = max(values(k))+1)))
#' plot(o, col = c(rgb(0,0,0,0), "red"), legend = FALSE, add = TRUE) # add origin
#'
#' a <- kissmigAccess(k, rel = TRUE)
#' plot(a*l, asp = 1.0, main = "Accessibility based on 'FOC', relative values",
#'      col = c("lightgrey", map.pal("viridis", n = max(values(k))+1)))
#' plot(o, col = c(rgb(0,0,0,0), "red"), legend = FALSE, add = TRUE) # add origin
#'
#' k <- kissmig(o, s, it = 150, type = "NOC")
#' plot(k*l, asp = 1.0, main = "Number of iteration steps with occurrences (NOC)",
#'      col = c("lightgrey", map.pal("viridis", n = max(values(k))+1)))
#' plot(o, col = c(rgb(0,0,0,0), "red"), legend = FALSE, add = TRUE) # add origin
#' }
#' @importFrom Rcpp evalCpp
#' @importFrom methods is
#' @import terra
#' @useDynLib kissmig, .registration = TRUE
#' @export kissmig

kissmig <- function(O, S = NULL, it, type = "FOC", signed = FALSE, pext = 1.0, pcor = 0.2, seed = NULL, n_threads = 1, n_random = 10000) {
  # check class of origin 'O' and read data
  if (!is(O, "SpatRaster")) stop("origin 'O' must be a SpatRaster of the terra package")
  if (dim(O)[3]!=1) stop("origin 'O' must be a SpatRaster with a single layer")
  ans <- O
  ov <- values(O)

  # check 'type'
  type <- toupper(type)
  ifelse(type %in% c("DIS", "FOC", "LOC", "NOC"),
    ty <- which(c("DIS", "FOC", "LOC", "NOC") == toupper(type)),
    stop("'type' must be 'DIS', 'FOC', 'LOC', or 'NOC'", call. = FALSE)
  )

  # check 'S', prepare suitability vector and dimensions
  ifelse(is.null(S),
    {
      sv <- rep(1.0, ncell(O))
      dh <- dim(O)
      warning("no suitability data found - globally set to 1.0", call. = FALSE)
    },
    {
      if (!(is(S, "SpatRaster"))) {
        stop("suitability 'S' must be a SpatRaster")
      }
      compareGeom(O, S) # stops if not the same extent, projection etc.
      sv <- as.vector(values(S)) # converts matrix in case of brick or stack
      dh <- dim(S)
    }
  )

  # check dimensions for cpp values
  # all dimensions have to be <= uint32_t MAX from c++
  # one layer has to be <= uint32_t MAX from c++
  # x < uint32Max
  # y < uint32Max
  # z < uint32Max
  # x * y < uint32Max
  uint32Max <- (2^32) - 1
  if (dh[1] > uint32Max) {
    stop(paste("the x dimension of a layer has to be <= (2^32)-1, found x =", dh[1]))
  } else if (dh[2] > uint32Max) {
    stop(paste("the y dimension of a layer has to be <= (2^32)-1, found y =", dh[2]))
  } else if (dh[3] > uint32Max) {
    stop(paste("the amount of layers in 'S' has to be <= (2^32)-1, found layers =", dh[3]))
  } else if ((dh[1] * dh[2]) > uint32Max) {
    stop(paste("the layer size has to be <= (2^32)-1, found x * y =", (dh[1] * dh[2]), ", for x =", dh[1], ", y =", dh[2]))
  }

  # check for too small dimensions on layer
  if (dh[1] < 3 || dh[2] < 3) {
    stop(paste("the smallest dimension per layer has to be at least >= 3, found x =", dh[1], ", y =", dh[2]))
  }

  # prevent reuse of the same random numbers per cell within or over a multiple of rows of one iteration
  # this is the case if columns * 11 are a multiple or a divider of n_random
  if (
    (((dh[2] * 11) %% n_random) == 0) ||
    ((n_random %% (dh[2] * 11)) == 0)
  ) {
    stop(paste0("'n_random' (", n_random, ") cannot be a multiple or a divisor of the amount of columns (", dh[2], ") multiplied by 11, this would generate undesirable pattern in the results due to internal random number usage\nPleas increase 'n_random'."))
  }

  # check for too large or small random numbers
  if (n_random < 10000) {
    warning(paste("amount of random numbers are smaller than the recommendation of 10000, found n_random =", n_random), call. = FALSE)
  } else if (n_random > uint32Max) {
    stop(paste("amount of random numbers have to be < (2^32)-1, found n_random =", n_random))
  }

  # check for too many iterations or layers
  # "(t * it) + l + 1" (required for "FOC","LOC") in (cpp) cannot exceed max value of uint32_t
  # requested iterations * number of layers + requested iterations (
  # +1 can be left since l will never reach it due to index by 0)
  if (
    !is.null(S) &&
      (((it * dh[3]) + it) > uint32Max)
  ) {
    stop(paste("the combination of ('it' * amount of layers in S + 'it') has to be <= (2^32)-1, found it = ", it, ", layers in S =", dh[3]))
  }

  # R RNG (seed isn't a parameter of kissmig_C)
  if (!is.null(seed)) set.seed(seed)

  # signed?
  si <- ifelse(signed, 1, 0)

  if (n_threads < 1) {
    n_threads <- 1
    warning("invalid thread count found, set to 1", call. = FALSE)
  }

  v <- kissmig_c(ov, sv, dh, it, pext, pcor, ty, si, n_threads, n_random)

  values(ans) <- v
  return(ans)
}
