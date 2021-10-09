#' @useDynLib driftSim
#' @import stats
#' @import utils
#' @import grDevices
#' @import graphics
#' @import ggplot2
#' @import Rcpp
.onUnload <- function(libpath){
    library.dynam.unload("driftSim",libpath)
}
