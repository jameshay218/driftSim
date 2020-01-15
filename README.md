# driftSim
> Short-sighted evolution of influenza cellular receptor binding avidity shapes influenza epidemic dynamics


[![Project Status: WIP - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)

This git repository contains all of the code needed to run the binding avidity adaptation project either through command line or the accompanying shiny app.

If you are using this app, you will likely be familiar enough with the code to get going, or have sufficient knowledge to work from the detailed vignettes listed here:

1. Summary of within-host aspect of the model - [the relationship between binding avidity, host immunity and transmission](https://jameshay218.github.io/driftSim/inst/doc/science.html).
2. Summary of [transmission model](https://jameshay218.github.io/driftSim/inst/doc/transmission_model.html).
3. Summary of [model output files](https://jameshay218.github.io/driftSim/inst/doc/outputs.html).
4. Technical [summary](https://jameshay218.github.io/driftSim/inst/doc/technical.html) of `driftSim`.

However, if you just want to get the app installed and running, see the minimal examples below.

## To do
1. Profile C++ code
2. Add reseeding code. This will take a sample of extant viruses at a given time, set this to one side and let it mutate, and then re-seed the epidemic with this virus at a later time
3. Virus phylogeny thinning. Every *n* steps, remove a number of viruses from the simulation, replacing each removed virus with its closests relative. This will slightly reduce accuracy, but should allow longer, bigger simulations to be run much quicker.
4. Change time step size in $\tau$-leap algorithm to be specified by the user. Currently set to 1 day time steps.
5. Full tests.

## Background
This project is a continuation of previously [published work](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3678328/). Here, we are interested in the hypothesis that within-host cellular binding avidity changes of influenza viruses play an important role in the evolutionary dynamics of influenza. Furthermore, we consider how antigenic drift interacts with binding avidity on an evolutionary basis. This package contains an implementation of an individual based model that incorporates binding avidity as a virus property, human host immunity, antigenic drift, and SIR dynamics in a host population.

## Installation
Installation is fairly straightforward. The package uses compiled C++ code, so you should have Rtools installed.

```r
devtools::install_github("jameshay218/driftSim")
```

## Usage 1 - command line
The primary method of using the simulation is through an R console. There are a number of inputs for the simulation. These have default values, but the user should make sure to explicitly specify them to ensure understanding of the underlying dynamics.
```r
library(driftSim)

## Example parameter values attached to the package:
## 1. flags (flags for simulation output saving)
## 2. hostpars - vector of parameters related to host population. See ?exampleParameters
## 3. viruspars - vector of parameters related to virus population. See ?exampleParameters
## 4. deltaVMat - matrix of gradients describing how binding avidity changes within a 
##    particular host in a single day. See ?calculate_deltaV_matrix
attach(exampleParameters)

## NOTE - the default parameters include random antigenic drift.
## To recover the model in the accompanying paper (with no drift),
## set viruspars["probMut"] to 0
viruspars["probMut"] <- 0

sim_duration <- 365 ## Duration of simulation in days
version <- 1
scenario_descriptions(1) ## Which version of the simulation to run (1-4)
output <- run_simulation(flags=flags, hostpars=hostpars, viruspars=viruspars, deltaVMat=deltaVMat,
                         iniKs=NULL,start=0,end=365,input_files="",output_files=c("SIR.csv","","","","",""),
                         VERBOSE=TRUE,scenario=version)
sir <- read.table("SIR.csv",header=FALSE,sep=",")
plot_SIR(sir)
```

## Usage 2 - shiny app
The shiny app is intended for the following three uses:
* Parameter exploration for binding avidity and within-host survival dynamics
* Interface to the Cpp simulation, allowing the user to vary the parameter inputs and view SIR dynamics from the 1-4 scenarios.
* View of the virus phylogenies using results from the simulation run. This will use Sean’s matlab code or a .exe
```r
library(driftSim)
runSimulationApp()
```

## License

GPL-3 © [James Hay &lt;james.hay13@imperial.ac.uk&gt;](https://github.com/).

