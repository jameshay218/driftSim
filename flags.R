## Host parameters to use
hostpars <- c(
  S0 = 299900, ## Starting susceptible population size
  I0 = 100, ## Initial infecteds
  R0 = 0, ## Initial recovered
  c = 0.7, ## Contact rate per day (ie. each person makes on average 0.7 contacts per day)
  mu = 1/(70*365), ## Birth and death rate per day
  w = 0.04, ## Immunity waning rate
  g = 1/3.3, ## Recovery rate
  iniBind = 0.45, ## Initial binding avidity of seed virus
  meanBoost = 8, ## Average (of poisson) antibody boosting upon recovery
  iniDist = 3, ## Initial antigenic distance of seed virus to "root"
  saveFreq = 5, ## How often the distribution of host immunity values are saved
  maxTitre = 25 ## Maximum attainable titre for hosts
  )

viruspars <- c(
  p=4, ## parameter to control degree by which changes in binding avidity affect probability of escape from immune response
  r=70, ## parameter to control degree by which previous exposure reduce probability of immune escape
  q=1 ,## parameter to control the shape of the relationship between binding avidity and immune escape (shift on the x-axis)
  a=0.7, ## controls rate of changes of relationship between probability of successful replication and change in binding avidity
  b=3, ## parameter to control the shape of the relationship between probability of successful replication and changes in binding avidity
  n=4, ## number of offspring per virus replication event
  v=1, ## number of virions initially transmitted
  probMut=0.1, ## Probability of an mutation in antigenic properties in a given infection
  expDist=1.3, ## Average of the exponential distribution from which the size of antigenic mutatations is drawn
  kc=0, ## Rate of binding avidity adaptation
  VtoD=0 ## Proportion of mediating how much antigenic change occurs given a unit of binding avidity change
)

flags <- c(
  SIR_flag=TRUE, ## Flag to save SIR dynamics
  voutput_phylogenetic_flag=FALSE, ## Flag to save virus information for phylogenetic tree analysis (ie. virus parent information)
  voutput_pairwise_flag=FALSE, ## Flag to save pairwise distance matrix
  time_flag=FALSE, ## Flag to record time taken for simulation
  save_state=FALSE, ## Flag to save the final state of the simulation
  input_flag_generated=FALSE, ## Flag to use specified file as input for simulation
  input_flag_saved=FALSE, ## Flag to use specified file as input for simulation
  save_k=FALSE ## Flag to save the distribution of hosts in each K value at regular intervals
)

#devtools::load_all("~/Documents/Binding Avidity/driftSim/")
deltaVMat <- calculate_deltaV_matrix(viruspars, 3, hostpars["maxTitre"], 1)

exampleParameters <- list(hostpars=hostpars,viruspars=viruspars,flags=flags,deltaVMat=deltaVMat)

output <- run_simulation(flags=flags, hostpars=hostpars, viruspars=viruspars,deltaVMat=deltaVMat,
                         iniKs=NULL,start=0,end=365,input_files="",output_files=c("SIR.csv","","","","",""),
                         VERBOSE=TRUE,scenario=3)
sir <- read.table("SIR.csv",header=FALSE,sep=",")
plot_SIR(sir)
