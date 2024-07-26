#' Run a simple species migration model
#'
#' \command{kissmig} runs a simple, raster-based, stochastic migration model to simulate species migration and
#' range shifts. It uses a geographic area of origin along with suitability maps to iteratively run a simple
#' 3x3 cell algorithm. Specifically, it allows for generating accessibility maps for easy integration of limited
#' migration in species distribution models (Nobis and Normand 2014, <\doi{10.1111/ecog.00930}>).
#' @usage kissmig(O, S=NULL, it, type='FOC', signed=FALSE, pext=1.0, pcor=0.2, seed=NULL)
#' @param O a single RasterLayer of the geographic origin
#' @param S a Raster* object of suitability, i.e., a RasterLayer, RasterStack, or RasterBrick
#' @param it number of iteration steps
#' @param type type of result: final distribution ('DIS'), iteration step of first occurrence
#' ('FOC'), iteration step of last occurrence ('LOC'), or number of iteration steps with occurrence ('NOC')
#' @param signed if TRUE, the sign indicates whether the cells was colonized (positive) or uncolonized (negative)
#' after the last iteration step
#' @param pext propability [0,1] a colonized cell becomes uncolonized between iteration steps, i.e., the species gets locally extinct
#' @param pcor propability [0,1] corner cells are considered in the 3x3 cell neighborhood
#' @param seed integer used to set the seed of the random number generator
#' @details
#' Starting from origin "O" \command{kissmig} simulates migration for "it" iteration steps in a heterogeneous environment
#' characterised by the suitability layer(s) "S". The colonized cells of the origin "O" have value 1, uncolonized
#' cells value 0. In case "S" consists of several suitability layers to cover environmental change, "it" is applied to each
#' layer. Suitability ranges between 0 (unsuitable) and 1 (suitability maximum). \command{kissmig} uses a 3x3 algorithm
#' for species spread/migration. All cells get exstinct before an iteration step with probability "pext", and for
#' a recolonization or new colonization event corner cells within the 3x3 neighborhood are considers
#' with probability "pcor" ("pcor"=0.2 produces more realistic circular spread patterns - see Nobis & Normand 2014).
#' For runtime optimization, signed results are generate for "signed"=TRUE, i.e, in addtion to
#' the result type 'FOC, 'LCO', or 'NOC', the sign indicates the final distribution ('DIS') with positive values
#' beeing colonized and negative values beeing previously colonized but uncolonized after the last iteration step.
#' To get reproducible results, the seed of the R random number generator can be set using the "seed" parameter.
#' @references \itemize{
#' \item{Nobis MP and Normand S (2014) KISSMig - a simple model for R to account for
#' limited migration in analyses of species distributions. \cite{Ecography} 37: 1282-1287. <\doi{10.1111/ecog.00930}>}
#' \item{KISSMig homepage <\url{https://purl.oclc.org/wsl/kissmig}>}}
#' @seealso \code{\link{kissmigAccess}, \link{kissmigOrigin}}
#' @examples
#' library(kissmig)
#'
#' # create a suitability map and geographic origin
#'
#' s <- kissmigDummyS(mean=12, sd=3)
#' o <- kissmigOrigin(s, x=8, y=44.5, size=0.5)
#' l <- s>=0 # land mask used for plotting kissmig results below
#' plot(s, asp=1.0, main='suitability + origin (in black)')
#' plot(o, col=c(NA,"black"), legend=FALSE, add=TRUE) # add origin
#'
#' # run kissmig with different type of output
#'
#' k <- kissmig(o, s, it=150, type='FOC')
#' plot(k*l, asp=1.0, main='First iteration step of occurrence (type="FOC")')
#' plot(o, col=c(NA,"black"), legend=FALSE, add=TRUE) # add origin
#' \donttest{
#' a <- kissmigAccess(k)
#' plot(a*l, asp=1.0, main='Accessibility based on "FOC", absolute values')
#' plot(o, col=c(NA,"black"), legend=FALSE, add=TRUE) # add origin
#'
#' a <- kissmigAccess(k, rel=TRUE)
#' plot(a*l, asp=1.0, main='Accessibility based on "FOC", relative values')
#' plot(o, col=c(NA,"black"), legend=FALSE, add=TRUE) # add origin
#'
#' k <- kissmig(o, s, it=150, type='DIS')
#' plot(k*l, asp=1.0, main='Final distribution (type="DIS")')
#' plot(o, col=c(NA,"black"), legend=FALSE, add=TRUE) # add origin
#'
#' k <- kissmig(o, s, it=150, type='LOC')
#' plot(k*l, asp=1.0, main='Last iteration step of occurrence (type="LOC")')
#' plot(o, col=c(NA,"black"), legend=FALSE, add=TRUE) # add origin
#'
#' k <- kissmig(o, s, it=150, type='NOC')
#' plot(k*l, asp=1.0, main='Number of iteration steps with occurrences (type="NOC")')
#' plot(o, col=c(NA,"black"), legend=FALSE, add=TRUE) # add origin
#' }
#'
#' @importFrom Rcpp evalCpp
#' @importFrom methods is
#' @import raster
#' @useDynLib kissmig, .registration = TRUE
#' @export kissmig

kissmig <- function(O, S=NULL, it, type='FOC', signed=FALSE, pext=1.0, pcor=0.2, seed=NULL) {

  # check class of origin 'O' and read data
  if (!is(O, 'RasterLayer')) stop("origin 'O' must be a RasterLayer")
  ans <- O
  ov  <- values(O)

  # check 'type'
  type <- toupper(type)
  ifelse(type %in% c('DIS','FOC','LOC','NOC'),
    ty <- which(c('DIS','FOC','LOC','NOC')==toupper(type)),
    stop("'type' must be 'DIS', 'FOC', 'LOC', or 'NOC'", call. = FALSE)
  )

  # check 'S', prepare suitability vector and dimensions
  ifelse(is.null(S),
    { sv <- rep(1.0, ncell(O))
      dh <- dim(O)
      warning('no suitability data found - globally set to 1.0', call. = FALSE)
    },{
      if (!(is(S, 'RasterLayer') | is(S, 'RasterStack') | is(S, 'RasterBrick'))) {
      stop("suitability 'S' must be a RasterLayer, RasterStack, or RasterBrick")
      }
      compareRaster(O,S)           # stops if not the same extent, projection etc.
      sv <- as.vector(values(S))   # converts matrix in case of brick or stack
      dh <- dim(S)
    })

  # R RNG (seed isn't a parameter of kissmig_C)
  if (!is.null(seed)) set.seed(seed)

  # signed?
  si <- ifelse(signed, 1, 0)

  v <- .Call('kissmig_c',
             as.double(ov),
             as.double(sv),
             as.integer(dh),
             as.integer(it),
             as.double(pext),
             as.double(pcor),
             as.integer(ty),
             as.integer(si),
			 PACKAGE='kissmig'
  )

  values(ans) <- v
  return(ans)

}
