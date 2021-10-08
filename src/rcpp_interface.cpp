//rcpp_functions.cpp
//[[Rcpp::plugins(cpp11)]]

#include <iostream>
#include <vector>
#include <time.h>
#include <memory>
#include <fstream>
#include <string>
#include <ctime>
#include <Rmath.h>
#include <Rcpp.h>
#include <vector>

#include "host.hpp"
#include "hostpopulation.hpp"

using namespace std;

//' Interface to the entire C++ driftSim simulation
//'
//' @param flags set of boolean flags for simulation output. See \code{\link{run_simulation}}
//' @param hostPopn vector of parameters related to the host population. See \code{\link{run_simulation}}
//' @param vlPars Matrix of parameters related to the iviral load kinetics. See \code{\link{run_simulation}}
//' @param day integer of simulation starting day
//' @param final_day integer of simulation final day (steps of 1 day)
//' @param output_files vector of output file names, as in \code{\link{run_simulation}}
//' @param VERBOSE boolean to indicate how much console output the simulation should have
//' @param scenario integer 1 to 4, indicating which version of the model is being run. See \code{\link{scenario_descriptions}} for a description of each scenario
//' @param callback pointer to enable call backs to R. Pretty much just for progress bar in the shiny app
//' @return the final virus ID from the simulation
//' @export
//' @useDynLib driftSim
//[[Rcpp::export]]
int run_simulation_cpp(Rcpp::IntegerVector flags, 
		       Rcpp::NumericVector hostPopn, 
		       Rcpp::NumericMatrix vlPars,  
		       Rcpp::NumericMatrix infectiousnessPars,
		       Rcpp::NumericMatrix crossImmunity,
		       int day, 
		       int final_day,
		       std::vector<std::string> output_files, 
		       bool VERBOSE,
		       SEXP callback)
{
  HostPopulation* hpop;

  // Set up call back options.
  Rcpp::Environment myEnv("package:driftSim");

  // Extract option flags
  bool save_SIR = flags[0];
  bool save_viruses = flags[1];
  bool use_time = flags[2];
  bool save_hosts = flags[3];

  // Extract host population parameters
  double S0 = hostPopn[0];
  double I0 = hostPopn[1];
  double R0 = hostPopn[2];
  double start_day = day;
  double contactRate = hostPopn[3];
  double mu = hostPopn[4];
  double wane = hostPopn[5];
  
  // Virus parameters
  double true_0 = hostPopn[6];
  double seed_variant = hostPopn[7];
  int n_variants = infectiousnessPars.nrow();
  //  Timing parameters
  int startT = clock();
  int endT = 0;
  
  double new_infected;

  // Set up streams for simulation output
  string filename;
  ofstream output (output_files[0]);
  ofstream voutput;
  ofstream houtput;

  // Call back function? Pretty much just for progress bar in the shiny app
  const bool has_callback = callback != R_NilValue;
  if (has_callback) {
    callback = Rcpp::as<Rcpp::Function>(callback);
  }

  // Set the virus parameters
  Virus::set_true_0(true_0);
  Virus::set_variants(n_variants);
  
  // Print all starting conditions
  if(VERBOSE){
    Rcpp::Rcout << "###################" << endl;
    Rcpp:: Rcout << "Starting simulation..." << endl << endl;
    Rcpp::Rcout << "Parameters used: " << endl << endl;
    Rcpp::Rcout << "S0: " << S0 << endl;
    Rcpp::Rcout << "I0: " << I0 << endl;
    Rcpp::Rcout << "R0: " << R0 << endl;
    Rcpp::Rcout << "contact: " << contactRate << endl;
    Rcpp::Rcout << "mu: " << mu << endl;
    Rcpp::Rcout << "wane: " << wane << endl;
  }

   // Generate a host population with these properties
   Rcpp::Rcout << "Simulating population" << endl;
   hpop = new HostPopulation(S0,I0, R0, start_day,contactRate, mu, wane, seed_variant);


  if(VERBOSE){
    Rcpp::Rcout << "Starting conditions: " << endl;
    hpop->printStatus();
    Rcpp::Rcout << endl;
  }

  /* ================= MAIN LOOP =============== */
  while(day <= final_day){
    // Move the simulation forward by one day
    new_infected = hpop->stepForward(day);
      
    if(VERBOSE){
      hpop->printStatus();
      Rcpp::Rcout << endl;
    }
    
    // Write SIR values to output file
    if(save_SIR) output << hpop->countSusceptibles() << "," << hpop->countInfecteds() << "," << hpop->countRecovereds() << "," << new_infected << endl;

    // Update the progress bar
    if(has_callback) {
      Rcpp::IntegerVector tmp = Rcpp::IntegerVector::create(
							    Rcpp::_["day"]=hpop->getDay(),
							    Rcpp::_["susceptibles"]=hpop->countSusceptibles(),
							    Rcpp::_["infecteds"]=hpop->countInfecteds(),
							    Rcpp::_["recovereds"]=hpop->countRecovereds());
      Rcpp::as<Rcpp::Function>(callback)(tmp);
    }
    day++;
  }
  /* =============================== */

  if(VERBOSE) Rcpp::Rcout << "Simulation finished" << endl;
  // Save the virus population characteristics to file
  if(save_viruses) hpop->writeViruses(voutput, output_files[1]);
  // Save the human host population characteristics to file
  if(save_hosts){
    hpop->writeHosts(houtput, output_files[2]);
    Rcpp::Rcout << endl;
  }

  delete hpop;

  // Close the SIR output stream
  if(save_SIR) output.close();
  
  if(use_time){
    endT = clock();
    Rcpp::Rcout << "Time elapsed: " << (endT-startT)/double(CLOCKS_PER_SEC) << " Seconds" << endl;
  }
  return Virus::getIDgenerator();
}

