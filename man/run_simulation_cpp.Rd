% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{run_simulation_cpp}
\alias{run_simulation_cpp}
\title{Interface to the entire C++ driftSim simulation}
\usage{
run_simulation_cpp(flags, hostPopn, virusPars, deltaVMat, iniKs, day,
  final_day, input_files, output_files, VERBOSE, scenario, callback)
}
\arguments{
\item{flags}{set of boolean flags for simulation output. See \code{\link{run_simulation}}}

\item{hostPopn}{vector of parameters related to the host population. See \code{\link{run_simulation}}}

\item{virusPars}{vector of parameters related to the influenza virions. See \code{\link{run_simulation}}}

\item{deltaVMat}{matrix of gradients for use in simulating adaptive binding avidity.}

\item{iniKs}{data frame of initial K values}

\item{day}{integer of simulation starting day}

\item{final_day}{integer of simulation final day (steps of 1 day)}

\item{input_files}{vector of input file names, as in \code{\link{run_simulation}}}

\item{output_files}{vector of output file names, as in \code{\link{run_simulation}}}

\item{VERBOSE}{boolean to indicate how much console output the simulation should have}

\item{scenario}{integer 1 to 4, indicating which version of the model is being run. See \code{\link{scenario_descriptions}} for a description of each scenario}

\item{callback}{pointer to enable call backs to R. Pretty much just for progress bar in the shiny app}
}
\value{
the final virus ID from the simulation
}
\description{
Interface to the entire C++ driftSim simulation
}
