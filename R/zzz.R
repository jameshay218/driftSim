#' @useDynLib driftSim
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import utils
#' @import grDevices
#' @import graphics
#' @import ggplot2
.onUnload <- function(libpath){
    library.dynam.unload("driftSim",libpath)
}
