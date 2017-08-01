## ----setup, include=FALSE------------------------------------------------
library(driftSim)
knitr::opts_chunk$set(echo = TRUE, fig.width=7,fig.height=6,
                      message=FALSE, results=FALSE,warning=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  ## Generate a uniform distribution of 1000 k values from 0 to 9
#  iniKs <- floor(runif(1000,0, 10))
#  
#  ## How many individuals do we want in the new simulation?
#  N <- 10000
#  
#  ## From the same distribution of iniKs, draw N k values
#  startingKs <- generateHostKDist_2(iniKs,N)

## ----eval=FALSE----------------------------------------------------------
#  hostK_file <- read.csv("hosts.csv")
#  iniKs <- hostK_file$hostK
#  flags["input_flag_generated"] <- TRUE
#  output <- run_simulation(flags=flags,
#                           hostpars=hostpars, viruspars=viruspars, deltaVMat=deltaVMat,
#                           iniKs=iniKs,
#                           start=0,end=365,input_files="",output_files=c("SIR.csv","","","","",""),
#                           VERBOSE=TRUE,scenario=version)

