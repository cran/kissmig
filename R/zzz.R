.onAttach <- function(libname, pkgname) {
  msg <- paste(
    paste("kissmig", utils::packageVersion("kissmig")),
    "- new in v2.0: speedup through parallel processing")
  packageStartupMessage(msg)
}
