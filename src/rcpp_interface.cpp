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


//' Solve viral kinetics
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector solve_viral_kinetics_cpp(Rcpp::NumericVector times, double tg, double tp, double to, double tw, double alpha,
                                             double infectiousnessMax, double infectiousnessGradient, double infectiousnessInflection){
    Virus* V = new Virus(times[0],tg, tp, to, tw, alpha,infectiousnessMax,infectiousnessGradient,infectiousnessInflection);
    Rcpp::NumericVector vls(times.size());
    
    for(int t = 0; t < times.size(); t++){
        vls[t] = V->calculateViralLoad(times[t]);
    }
    return vls;
}
//' Solve infectiousness
//'
//' @export
//[[Rcpp::export]]
Rcpp::NumericVector solve_infectiousness_cpp(Rcpp::NumericVector times, double tg, double tp, double to, double tw, double alpha,
                                             double infectiousnessMax, double infectiousnessGradient, double infectiousnessInflection){
    Virus* V = new Virus(times[0],tg, tp, to, tw, alpha,infectiousnessMax,infectiousnessGradient,infectiousnessInflection);
    Rcpp::NumericVector infectiousness(times.size());
    
    for(int t = 0; t < times.size(); t++){
        infectiousness[t] = V->getInfectiousness(times[t]);
    }
    return infectiousness;
}




//' Interface to the entire C++ driftSim simulation
//'
//' @param flags set of boolean flags for simulation output. See \code{\link{run_simulation}}
//' @param hostPopn vector of parameters related to the host population. See \code{\link{run_simulation}}
//' @param seeds Matrix of parameters related to epidemic seeding. See \code{\link{run_simulation}}
//' @param vlPars Matrix of parameters related to the viral load kinetics. See \code{\link{run_simulation}}
//' @param infectiousnessPars Matrix of parameters related to the variant-specific infectiousness. See \code{\link{run_simulation}}
//' @param crossImmunity Matrix of parameters related to variant cross immunity. See \code{\link{run_simulation}}
//' @param day double of simulation starting day
//' @param final_day double of simulation final day (steps of 1 day)
//' @param tstep double of simulation time steps
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
		       Rcpp::NumericMatrix seeds,
		       Rcpp::NumericMatrix vlPars,  
		       Rcpp::NumericMatrix infectiousnessPars,
		       Rcpp::NumericMatrix crossImmunity,
		       Rcpp::NumericVector ageDistribution,
		       Rcpp::NumericVector infectivityCurve,
		       Rcpp::NumericVector susceptibilityCurve,
		       double day, 
		       double final_day,
		       double tstep,
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
  double R0 = hostPopn[1];
  double start_day = day;
  double contactRate = hostPopn[2];
  double mu = hostPopn[3];
  double wane = hostPopn[4];
  int n_tests = hostPopn[6];
  
  // First seeding
  double seed_time = seeds(0,0);
  int seed_variant = seeds(0,1);
  int I0 = seeds(0,2);
  
  // Virus parameters
  double true_0 = hostPopn[5];
  
  int n_variants = infectiousnessPars.nrow();
  //  Timing parameters
  int startT = clock();
  int endT = 0;
  
  double new_infected;

  // Set up streams for simulation output
  string filename;
  ofstream output (output_files[0]);
  ofstream output_tests (output_files[1]);
  ofstream voutput;
  ofstream houtput;

  // Call back function? Pretty much just for progress bar in the shiny app
  const bool has_callback = callback != R_NilValue;
  if (has_callback) {
    callback = Rcpp::as<Rcpp::Function>(callback);
  }

  // Set the virus parameters
  Virus::set_true_0(true_0);
  Virus::set_generator(0);
  Virus::set_variants(n_variants);
  Virus::initiateVariantGenerator();
  Virus::set_vlPars(vlPars);
  Virus::set_infectiousnessPars(infectiousnessPars);
  Virus::set_crossImmunity(crossImmunity);
  
  // Print all starting conditions
  if(VERBOSE){
    Rcpp::Rcout << "###################" << endl;
    Rcpp::Rcout << "Starting simulation..." << endl << endl;
    Rcpp::Rcout << "Parameters used: " << endl << endl;
    Rcpp::Rcout << "S0: " << S0 << endl;
    Rcpp::Rcout << "I0: " << I0 << endl;
    Rcpp::Rcout << "R0: " << R0 << endl;
    Rcpp::Rcout << "contact: " << contactRate << endl;
    Rcpp::Rcout << "mu: " << mu << endl;
    Rcpp::Rcout << "wane: " << wane << endl;
    Rcpp::Rcout << "Seed variant: " << seed_variant << endl;
  }

   // Generate a host population with these properties
   Rcpp::Rcout << "Simulating population" << endl;
  
   hpop = new HostPopulation(S0,I0, R0, start_day,contactRate, mu, wane, seed_variant, 
                             infectivityCurve, susceptibilityCurve,ageDistribution);
   Rcpp::Rcout << "Host population created" << std::endl;

  if(VERBOSE){
    Rcpp::Rcout << "Starting conditions: " << endl;
    hpop->printStatus();
    Rcpp::Rcout << endl;
  }

  /* ================= MAIN LOOP =============== */
   
   int seed_index = 1;
  while(day <= final_day){
      // Check for next seeding
      if(seed_index < seeds.nrow()){
        seed_time = seeds(seed_index,0);
        seed_variant = seeds(seed_index,1);
        I0 = seeds(seed_index,2);
        
        if(day == seed_time){
            Rcpp::Rcout << "Seeding!" << std::endl;
            hpop->seed(seed_variant, I0);
        }
      }
      
    // Move the simulation forward by one day
    new_infected = hpop->stepForward(day);
      
    if(VERBOSE){
      hpop->printStatus();
      Rcpp::Rcout << endl;
    }
    
    // Write SIR values to output file
    if(save_SIR) output << hpop->countSusceptibles() << "," << hpop->countInfecteds() << "," << hpop->countRecovereds() << "," << new_infected << endl;

    // Sample and save testing round
    hpop->testingRound(output_tests, n_tests);
    
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
  if(save_viruses) hpop->writeViruses(voutput, output_files[2]);
  // Save the human host population characteristics to file
  if(save_hosts){
    hpop->writeHosts(houtput, output_files[3]);
    Rcpp::Rcout << endl;
  }

  delete hpop;

  // Close the SIR output stream
  if(save_SIR) output.close();
  
  output_tests.close();
  
  if(use_time){
    endT = clock();
    Rcpp::Rcout << "Time elapsed: " << (endT-startT)/double(CLOCKS_PER_SEC) << " Seconds" << endl;
  }
  return Virus::getIDgenerator();
}

