## ----setup, include=FALSE, eval=TRUE,echo=FALSE--------------------------
library(driftSim)
knitr::opts_chunk$set(echo = TRUE, fig.width=7,fig.height=6,
                      message=FALSE, results=FALSE,warning=FALSE)



## ---- eval=TRUE, message=FALSE,warnings=FALSE----------------------------
## Basic example of producing these host files
attach(exampleParameters)

sim_duration <- 365 ## Duration of simulation in days
version <- 4
flags <- as.logical(c(1,1,1,0,1,0,0,1))
scenario_descriptions(4) ## Which version of the simulation to run (1-4)
viruspars["kc"] <- 0.2
output <- run_simulation(flags=flags, 
                         hostpars=hostpars, 
                         viruspars=viruspars, 
                         deltaVMat=deltaVMat,
                         iniKs=NULL,
                         start=0,end=365,
                         input_files="",
                         output_files=c("SIR.csv","voutput_characteristics.csv","voutput_distances.csv","hosts.csv","hostKs.csv"),                                   VERBOSE=TRUE,scenario=4)

## ------------------------------------------------------------------------
sir <- read.table("SIR.csv",header=FALSE,sep=",")
plot_SIR(sir)

## ----fig.width=7,fig.height=10-------------------------------------------
## Read in virus output data
tmp <- data.table::fread("voutput_characteristics.csv",data.table=TRUE)
dt <- data.table::data.table(tmp)

## Get mean and 95% CI of distance to root over time
means <- data.frame(dt[,list(mean=mean(distRoot)),by=birth])
lower <- data.frame(dt[,list(lower=quantile(distRoot,c(0.025))),by=birth])
upper <- data.frame(dt[,list(upper=quantile(distRoot,c(0.975))),by=birth])

## Also take a random virus that was extant at the
## end of the simulation, and use this to look at
## an example antigenic change profile
virus_data <- data.table::fread("voutput_characteristics.csv",data.table=FALSE)
## Note that extant viruses are given a death time of -1
samps <- virus_data[(virus_data$death == -1 & virus_data$birth > 0),]
i <- sample(nrow(samps),1)
startVid <- virus_data[virus_data$vid==samps[i,"vid"],"vid"]
index <- 1
nextVid <- 1

## Trace back through this virus lineage until we reach the
## root virus.
x <- y <- z <- NULL
while(nextVid != 0){
  nextVid <- virus_data[virus_data$vid==startVid,"parentid"]
  x[index] <- virus_data[virus_data$vid==startVid,"birth"]
  y[index] <- virus_data[virus_data$vid==startVid,"bindingAvidityFinal"]
  z[index] <- virus_data[virus_data$vid==startVid,"distRoot"]
  index <- index + 1
  startVid <- nextVid
}
dat <- data.frame(time = x, V=y ,delta=z )
dat1 <- dat[with(dat,order(time,delta)),]
dat1 <- dat1[dat1$time > 0,]

par(mfrow=c(2,1))
plot(means,type='l',lwd=1,col="blue", ylim=c(0,40),
     xlab="Time (days)",ylab="Average Antigenic Distance to Root",
     main="Mean and 95% CI of antigenic distances to root over time")
lines(lower,col="blue",lwd=0.5)
lines(upper,col="blue",lwd=0.5)

plot(dat1$delta~dat1$time,type='l',col="red",
     xlab="Time (days)", ylab="Antigenic distance from root",
     main="Example antigenic distance over time of an extant lineage")

## ----fig.width=6,fig.height=6--------------------------------------------
voutput2 <- data.table::fread("voutput_distances.csv",data.table=FALSE)
## An uninformative plot, but just to show an example use
image(as.matrix(voutput2))

## ----fig.width=6,fig.height=4--------------------------------------------
hostKDist <- generateHostKDist("hosts.csv",100000)
hist(hostKDist)

## ----fig.width=7,fig.height=4--------------------------------------------
hostKs <- read.table("hostKs.csv",sep=",")
day <- 25*hostpars["saveFreq"]

barplot(as.numeric(hostKs[25,1:30]),
        xlab="k",ylab="Frequency",
        main=paste0("Frequency of host k values at day ",day))

