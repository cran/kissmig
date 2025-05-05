.onAttach <- function(libname, pkgname) {
  msg <- paste(
    paste("kissmig", utils::packageVersion("kissmig")),
    "- major new feature: speedup through parallel processing")
  packageStartupMessage(msg)
}
