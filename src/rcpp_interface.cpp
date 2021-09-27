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
//' @param virusPars vector of parameters related to the influenza virions. See \code{\link{run_simulation}}
//' @param deltaVMat matrix of gradients for use in simulating adaptive binding avidity.
//' @param iniKs data frame of initial K values
//' @param day integer of simulation starting day
//' @param final_day integer of simulation final day (steps of 1 day)
//' @param input_files vector of input file names, as in \code{\link{run_simulation}}
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
		       Rcpp::NumericVector virusPars,  
		       Rcpp::NumericMatrix deltaVMat, 
		       SEXP iniKs,
		       int day, 
		       int final_day,
		       std::vector<std::string> input_files, 
		       std::vector<std::string> output_files, 
		       bool VERBOSE,
		       int scenario,
		       SEXP callback)
{
  int MAXK = 301;
  HostPopulation* hpop;

  // Set up call back options.
  Rcpp::Environment myEnv("package:driftSim");
  Rcpp::Function generateHostKDist = myEnv["generateHostKDist"];
  Rcpp::Function writeCSV = myEnv["writeCSV"];
  Rcpp::NumericVector startingKs;

  // Extract option flags
  bool save_SIR = flags[0];
  bool save_viruses = flags[1];
  bool save_pairwise_viruses = flags[2];
  bool use_time = flags[3];
  bool save_hosts = flags[4];
  bool import_start_generate = flags[5];
  bool import_start_saved = flags[6];
  bool save_hostKs = flags[7];

  // Extract host population parameters
  double S0 = hostPopn[0];
  double I0 = hostPopn[1];
  double R0 = hostPopn[2];
  double start_day = day;
  double contactRate = hostPopn[3];
  double mu = hostPopn[4];
  double wane = hostPopn[5];
  double gamma = hostPopn[6];
  double iniBinding = hostPopn[7];
  double meanBoost = hostPopn[8];
  double iniDist = hostPopn[9];
  double kSaveFreq = hostPopn[10];
  double maxTitre = hostPopn[11];

  // Extract virus population parameters
  Rcpp::NumericMatrix allCounts(((final_day-day)/kSaveFreq)+1,MAXK);
  double p = virusPars[0];
  double r = virusPars[1];
  double q = virusPars[2];
  double a = virusPars[3];
  double b = virusPars[4];
  double n = virusPars[5];
  double v = virusPars[6];
  double probMut = virusPars[7]; // Probability of an antigenic mutation upon infection
  double expDist = 1/virusPars[8]; // Exponential distribution rate of antigenic mutation size
  double kc = virusPars[9];
  double VtoD = virusPars[10];

  //  Timing parameters
  int index = 0;
  int startT = clock();
  int endT = 0;
  
  double new_infected;

  // Set up streams for simulation output
  string filename;
  ofstream output (output_files[0]);
  ofstream voutput;
  ofstream voutput2;
  ofstream houtput;

  // Call back function? Pretty much just for progress bar in the shiny app
  const bool has_callback = callback != R_NilValue;
  if (has_callback) {
    callback = Rcpp::as<Rcpp::Function>(callback);
  }


  // Set the virus parameters
  Virus::set_p(p);
  Virus::set_r(r);
  Virus::set_q(q);
  Virus::set_a(a);
  Virus::set_b(b);
  Virus::set_n(n);
  Virus::set_v(v);
  Virus::set_prob_mut(probMut);
  Virus::set_exp_dist(expDist);
  Virus::set_kc(kc);
  Virus::set_VtoD(VtoD);
  Virus::set_scenario(scenario);
  Virus::set_generator(1);
  Virus::set_delta(iniDist);
  Virus::set_deltaVMat(deltaVMat);

  Host::changeMeanBoost(meanBoost);
  Host::set_maxTitre(maxTitre);

  Rcpp::Rcout << "Max titre: " << maxTitre << endl;
  Rcpp::Rcout << "Mutation rate: " << probMut << " and " << expDist << endl;
  
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
    Rcpp::Rcout << "gamma: " << gamma << endl;
    Rcpp::Rcout << "iniBinding: " << iniBinding << endl;
    Rcpp::Rcout << "Scenario: " << scenario << endl << endl;
  }

  // If a distribution of host K values has been given and the flag for using saved K values is specified
  // Generate a host population with these properties
  if(iniKs != R_NilValue && import_start_saved){
    Rcpp::Rcout << "Simulating population using provided immunity values" << endl;
    hpop = new HostPopulation(S0, I0,R0,start_day,contactRate,mu,wane,gamma,iniBinding,iniDist, Rcpp::as<Rcpp::NumericVector>(iniKs));	
  } else if(import_start_generate){
    // Otherwise, if the flag for generating a new host K distribution has been specified, 
    // call back to R to get this distribution
    Rcpp::Rcout << "Simulating population using provided immunity distribution" << endl;
    startingKs = callFunction(input_files[0],(S0 + I0 + R0), generateHostKDist);
    hpop = new HostPopulation(S0, I0,R0,start_day,contactRate,mu,wane,gamma,iniBinding,iniDist, startingKs);		      
  } else {
    // Otherwise, create an entirely new host population (ie. no immunity)
    Rcpp::Rcout << "Simulating population, no immunity" << endl;
    hpop = new HostPopulation(S0,I0, R0, start_day,contactRate, mu, wane, gamma,iniBinding, iniDist);
  }

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
    // If saving K values are being saved over time, record these
    if(save_hostKs && day%(int)kSaveFreq == 0){
      // Extract the distribution of host K values
      Rcpp::NumericVector tmp1 = hpop->getHostKDist();
      // Count the number in each K category
      Rcpp::NumericVector tmptmp = countKs(tmp1, MAXK);
      allCounts(index,Rcpp::_) = tmptmp;
      index++;
    }

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
  if(save_viruses) hpop->writeViruses(voutput, output_files[1], FALSE);
  // Save the virus pairwise antigenic distances to file
  if(save_pairwise_viruses) hpop->virusPairwiseMatrix(voutput2, output_files[2],1000);
  // Save the human host population characteristics to file
  if(save_hosts){
    hpop->writeHosts(houtput, output_files[3]);
    Rcpp::Rcout << endl;
  }

  // Finally write the hhost population K counts to file
  if(save_hostKs) writeCSV(allCounts,output_files[4]);

  delete hpop;

  // Close the SIR output stream
  if(save_SIR) output.close();
  
  if(use_time){
    endT = clock();
    Rcpp::Rcout << "Time elapsed: " << (endT-startT)/double(CLOCKS_PER_SEC) << " Seconds" << endl;
  }
  return Virus::getIDgenerator();
}

