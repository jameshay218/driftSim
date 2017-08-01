# driftSim
> Short-sighted evolution of influenza cellular receptor binding avidity shapes influenza epidemic dynamics


[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

This git repository contains all of the code needed to run the binding avidity adaptation project either through command line or the accompanying shiny app.

If you are using this app, you will likely be familiar enough with the code to get going already, or have sufficient knowledge to work from the detailed [vignette](). However, if you just want to get the app installed and running, see the minimal examples below.

The essential scientific documentation of the model and what it aims to achieve can be found [here]().



## Background
This project is a continuation of previously [published work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3678328/). Here, we are interested in the hypothesis that within-host cellular binding avidity changes of influenza viruses play an important role in the evolutionary dynamics of influenza. Furthermore, we consider how antigenic drift interacts with binding avidity on an evolutionary basis. This package contains an implementation of an inidividual based model that incorporates binding avidity as a virus property, human host immunity, antigenic drift, and SIR dynamics in a host population.

## Installation
Installation is fairly straightforward. The package uses compiled C++ code, so you should have Rtools installed.

```r
devtools::install_github("jameshay218/driftSim")
```

## Usage 1 - command line
The first way of using the simulation is through an R console. There are a number of inputs for the simulation. These have default values, but the user should make sure to explicitly specify each to ensure understanding of the underlying dynamics.
```r
library(driftSim)

## Example parameter values attached to the package:
## 1. flags (flags for simulation output saving)
## 2. hostpars - vector of parameters related to host population. See vignettes
## 3. viruspars - vector of parameters related to virus population. See vignettes
## 4. deltaVMat - matrix of gradients describing how binding avidity changes within a 
##    particular host in a single day
attach(exampleParameters)

sim_duration <- 365 ## Duration of simulation in days
version <- 1
scenario_descriptions(1) ## Which version of the simulation to run (1-4)
output <- run_simulation(flags=flags, hostpars=hostpars, viruspars=viruspars, deltaVMat=deltaVMat,
                         iniKs=NULL,start=0,end=365,input_files="",output_files=c("SIR.csv","","","","",""),
                         VERBOSE=TRUE,scenario=3)
sir <- read.table("SIR.csv",header=FALSE,sep=",")
plot_SIR(sir)
```

## Usage 2 - shiny app
The shiny app is intended for the following three uses:
* Parameter exploration for binding avidity and within-host survival dynamics
* Interface to the Cpp simulation, allowing the user to vary the parameter inputs and view SIR dynamics from the 3-4 scenarios.
* View of the virus phylogenies using results from the simulation run. This will use Sean’s matlab code or a .exe
```r
library(driftSim)
runSimulationApp()
```

## License

GPL-2 © [James Hay &lt;james.hay13@imperial.ac.uk&gt;](https://github.com/).

