#' Main simulation interface
#' 
#' @param flags Vector of boolean flags for simulation outputs. All of the output file names are specified in the output_files arugment.
#' 1) save_SIR: if TRUE, saves the SIR dynamics
#' 2) save_viruses: if TRUE, saves information on the simulation viruses. From this file, phylogenetic reconstruction should be possible.
#' 3) save_pairwise_viruses: if TRUE, saves further virus information, spefically a matrix of pairwise antigenic distances between all viruses
#' 4) use_time: if TRUE, records the duration of the simulatoin run
#' 5) save_hosts: if TRUE, saves the properties of the simulation host population.
#' 6) import_start_generate: if TRUE, generates starting immunity (K) values for the host population. Requires iniK to be valid/
#' 7) import_start_saved: if TRUE, uses the provided iniKs vector as the actual starting immunity profile of the host population.
#' 8) save_hostsKs: if TRUE, saves the distribution of host immunity over time to file
#' @param hostpars vector of parameters relating to the host population. In order, these are: S0, I0, R0, contact rate, birth/death rate, temporary immunity waning rate (duration of recovered period), infected recovery rate, initial binding avidity of all simulated viruses, mean level of antibody boosting following recovered, initial antigenic distance of virus population to base virus, how often to save the host population K distribution, and the maximum achievable antibody titre by a given host. SEE THE VIGNETTES or \code{\link{exampleParameters}} FOR MORE DETAIL.
#' @param viruspars vector of parameters related to the virus population. These are, in order: p, r, q, a, b, n, v, probability of an antigenic mutation upon infection, mean of the exponential distribution of mutation sizes, kc and VtoD. SEE VIGNETTES or \code{\link{exampleParameters}} FOR DETAILS
#' @param deltaVMat A two dimensional matrix giving binding avidity changes for different immunity/binding avidity levels
#' @param start Simulation start day
#' @param end Simulation end day
#' @param output_files Vector of output file names. Note that these should be csv files. In order: 1) location of the SIR dynamics; 2) location of virus characteristics output; 3) location of pairwise antigenic distances; 4) Where to save the entire host population characteristics; 5) where to save the host population K distribution over time
#' @param VERBOSE If TRUE, prints additional simulation output
#' @param scenario Which version of the simulation to use. See \code{\link{scenario_descriptions}}
#' @param callback Leave this - this is just used by the shiny app for the progress bar.
#' @return Returns the ID of the last generated virus
#' @export
#' @useDynLib driftSim
#' @examples
#' attach(exampleParameters)
#' sim_duration <- 365 ## Duration of simulation in days
#' version <- 1
#' scenario_descriptions(1) ## Which version of the simulation to run (1-4)
#' output <- run_simulation(flags=flags, hostpars=hostpars, viruspars=viruspars, deltaVMat=deltaVMat,
#'                       iniKs=NULL,start=0,end=365,input_files="",output_files=c("SIR.csv","","","","",""),
#'                       VERBOSE=TRUE,scenario=3)
#' sir <- read.table("SIR.csv",header=FALSE,sep=",")
run_simulation <- function(
                           flags=c(1,0,0,0,0,0,0,0),
                           hostpars=c(90000,100,0,1.5,1/(40*365),1/25,0.333,0.8,10,0,10, 10),
                           viruspars=c(3, 70, 1, 0.7, 3, 2, 2, 0.1, 1, 0.5, 1000),
                           deltaVMat = matrix(ncol=2,nrow=2),
                           iniKs = NULL,
                           start=0,
                           end=100,
                           input_files=c("hosts.csv"),
                           output_files = c("SIR.csv","voutput1.csv","voutput2.csv","hosts.csv", "hostKs.csv"),
                           VERBOSE=TRUE,
                           scenario=1,
                           callback=NULL){
    if(scenario != 1 && is.null(deltaVMat)){
        print("Error - attempting binding avidity change with no deltaV matrix. Aborting")
        return(0)
    }
    print(input_files)
    return(run_simulation_cpp(flags,hostpars,viruspars, deltaVMat, iniKs, start,end,input_files,output_files,VERBOSE,scenario,callback))
}


##' Generate distribution of host immunity
##'
##' Given a vector of host immune values (K), creates a frequency distribution and then draws N samples to give a new host immunity distribution.
##' @param hostKs vector of K values
##' @param N number of samples to generate
##' @return a vector of random K values of size N
##' @export
##' @seealso \code{\link{generateHostKDist}}
generateHostKDist_2 <- function(hostKs, N){
  countHostK <- count(hostKs)
  freqs <- countHostK$freq/sum(countHostK$freq)
  cumSumK <- cumsum(freqs)
  startingKs <- generateKSamples(cumSumK, N)
  startingKs <- startingKs + 1
  return(countHostK$x[startingKs])
}

##' Generate distribution of host immunity from file
##'
##' Given a file which contains a column called "hostK", uses this column to generate a host K distribution as in \code{\link{generateHostKDist_2}}
##' @param hostFile the file name to be read in (a csv file)
##' @param N number of samples to generate
##' @return a vector of random K values of size N
##' @export
##' @useDynLib driftSim
##' @seealso \code{\link{generateHostKDist_2}}
generateHostKDist <- function(hostFile, N){
    hostDat <- read.csv(hostFile)
    hostKs <- hostDat$hostK
    final <- generateHostKDist_2(hostKs,N)
    return(final)
}

##' Launches the shiny app interface
##'
##' @export
##' @useDynLib driftSim
runSimulationApp <- function(){
        appDir <- system.file("shiny-examples", "driftSimApp", package = "driftSim")
        if (appDir == "") {
            stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
        }
        shiny::runApp(appDir, display.mode = "normal")
}


##' @export
##' @import Rcpp
##' @useDynLib driftSim
writeCSV <- function(tab, filename){
    write.table(tab,filename,col.names=FALSE,row.names=FALSE,sep=",")
}

##' Descriptions of model scenarios
##'
##' Returns a text description of a given scenario 1-4
##' @param scenario an integer 1-4 for the desired scenario
##' @return the text description of the given scenario
##' @export
scenario_descriptions <- function(scenario){
    if(!(scenario %in% 1:4)) return("Scenario does not exist")
    descriptions <- c("Random anigenic drift only. Binding avidity fixed",
                      "No drift; adaptive binding avidity; adaptive antigenic change. Binding avidity adaptation and antigenic change proportional to binding avidity change",
                      "Random drift, binding avidity adaptation and antigenic change from binding avidity change",
                      "Random drift and binding avidity adaptation (no binding avidity related antigenic change)")
    return(descriptions[scenario])
}

#' SIR dynamics plot
#'
#' Plots SIR dynamics for a given matrix of SIR time series data
#' @param dat the two dimensional matrix of S, I and R data
#' @return a ggplot2 object plot
#' @export
plot_SIR <- function(dat){
    dat <- cbind(seq(1,nrow(dat),by=1),dat)
    colnames(dat) <- c("t","S","I","R")
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

#' Calculate binding avidity adaptation matrix
#'
#' Calculates the gradient of binding avidity adaptation for different levels of host immunity.
#' This is useful because we need to solve some ODEs for each combination of host immunity and
#' each possible binding avidity value ie. precomputation
#' @param viruspars the vector of virus pars as taken by \link{\code{run_simulation}}
#' @param maxV the maximum binding avidity to calculate for
#' @param maxK the maximum achievable antibody titre in which the virus will have to adapt
#' @param step over how much time does this adaptation take place? usually 1 day
#' @export
calculate_deltaV_matrix <- function(viruspars, maxV=3,maxK=25, step = 1){
    message("Calculating deltaV matrix...")

    p <- viruspars["p"]
    r <- viruspars["r"]
    b <- viruspars["b"]
    a <- viruspars["a"]
    kc <- viruspars["kc"]
    q <- viruspars["q"]
    
    ## Differential equation to solve
    difeq <- function(t, V, params){
        x <- V
        j <- params[1]
        
        immK <- r*j
        if(immK <0) immK <- 0

        ## This is calculating the rate of binding avidity adaptation assuming
        ## that the virus tends towards to peak of f(x) and g(x), where these
        ## are the probabilities of successful replication and
        ## successful transmission. d/dx (f(x)g(x)) = g'(x)f(x) + g(x)f'(x)
        f_x <- (1-exp(-p*(x+q)))^(immK)
        g_x <- exp(-a*(x^b))
        f_dx <- p*(immK)*((1-exp(-p*(V+q)))^(immK-1))*(exp(-p*(x+q)))
        g_dx <- -a*b*x^(b-1)*exp(-a*(x^b))
        
        dV <- f_x*g_dx + g_x*f_dx
        return(list(dV))
    }
    V = seq(0,maxV,by=0.01)
    immKs = seq(0,maxK,by=0.1)
    allV <- matrix(nrow=length(immKs),ncol=length(V))
    for(j in 1:length(immKs)){
        message(cat("Current j: ", immKs[j]))
        for(i in 1:length(V)){
            deltaV <- deSolve::ode(y=c(V[i]),seq(0,step,by=1/40),difeq,c(immKs[j],viruspars))
            allV[j,i] <- deltaV[length(deltaV)] - V[i]
        }
    }
    return(allV)
}
