#' Main simulation interface
#' 
#' @param flags Vector of boolean flags for simulation outputs. All of the output file names are specified in the output_files arugment.
#' 1) save_SIR: if TRUE, saves the SIR dynamics
#' 2) save_viruses: if TRUE, saves information on the simulation viruses.
#' 3) use_time: if TRUE, records the duration of the simulation run
#' 4) save_hosts: if TRUE, saves the properties of the simulation host population.
#' @param hostpars vector of parameters relating to the host population. In order, these are: S0, I0, R0, contact rate, birth/death rate, temporary immunity waning rate (duration of recovered period), true 0 viral load value, and the ID of the seed variant
#' @param vlPars pre-computed matrix of random viral load kinetics parameters
#' @param infectiousnessPars matrix of infectiousness parameters
#' @param crossImmunity matrix of cross immunity
#' @param start Simulation start day
#' @param end Simulation end day
#' @param output_files Vector of output file names. Note that these should be csv files. In order: 1) location of the SIR dynamics; 2) location of virus characteristics output; 3) Where to save the entire host population characteristics.
#' @param VERBOSE If TRUE, prints additional simulation output
#' @param callback Leave this - this is just used by the shiny app for the progress bar.
#' @return Returns the ID of the last generated virus
#' @export
#' @useDynLib driftSim
run_simulation <- function(
    flags=c(1,0,0,0),
    hostpars=c(90000,100,0,1.5,1/(40*365),1/25,0,0),
    vlPars,
    infectiousnessPars,
    crossImmunity,
    start=0,
    end=100,
    output_files = c("SIR.csv","voutput.csv","hosts.csv"),
    VERBOSE=TRUE,
    callback=NULL){
    return(run_simulation_cpp(flags,hostpars,vlPars,infectiousnessPars,crossImmunity,
                              start, end, output_files,VERBOSE,callback))
}


#' SIR dynamics plot
#'
#' Plots SIR dynamics for a given matrix of SIR time series data
#' @param dat the two dimensional matrix of S, I and R data
#' @return a ggplot2 object plot
#' @export
plot_SIR <- function(dat){
    dat <- cbind(seq(1,nrow(dat),by=1),dat)
    colnames(dat) <- c("t","S","I","R","incidence")
    dat <- reshape2::melt(dat,id="t")
    colnames(dat) <- c("t","Population","value")
    SIR_plot <- ggplot() +
        geom_line(data=dat,aes(x=t,y=value,colour=Population,group=Population)) +
        xlab("Time (days)") +
        ylab("Number Individuals") +
        scale_y_continuous(expand=c(0,0))+
        scale_x_continuous(expand=c(0,0))+
        theme(
            text=element_text(colour="gray20",size=14),
            plot.title=element_text(size=28),
            legend.text=element_text(size=20,colour="gray20"),
            legend.title=element_text(size=20,colour="gray20"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line=element_line(colour="gray20"),
            axis.line.x = element_line(colour = "gray20"),
            axis.line.y=element_line(colour="gray20"),
            axis.text.x=element_text(colour="gray20"),
            panel.background=element_blank(),
            axis.text.y=element_text(colour="gray20"))
    return(SIR_plot)
}