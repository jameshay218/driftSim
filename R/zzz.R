#' @useDynLib driftSim
#' @importFrom Rcpp evalCpp
#' @importFrom deSolve ode
#' @import stats
#' @import utils
#' @import grDevices
#' @import graphics
#' @import ggplot2
.onUnload <- function(libpath){
    library.dynam.unload("driftSim",libpath)
}
