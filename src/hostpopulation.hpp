#ifndef HOSTPOPULATION_HPP
#define HOSTPOPULATION_HPP

#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <Rcpp.h>
#include <map>

class Host;
class Virus;

class HostPopulation{
private:
  std::vector<Host*> susceptibles;
  std::vector<Host*> infecteds;
  std::vector<Host*> recovereds;
  std::vector<Host*> dead;

  std::vector<Host*> new_infecteds;
  std::vector<Host*> new_recovereds;
  std::vector<Host*> new_susceptibles;
  std::vector<Host*> new_births;
  
  std::vector<Virus*> seed_viruses;

  double day;
  double contactRate; // How many contacts made on average per time step?
  double mu; // Birth and death rate
  double wane; // Rate of losing full immunity
  
  double inf_wane_rate;
  double vacc_wane_rate;
  double max_vacc_immunity;
  
  double vacc_prob;
  
  Rcpp::NumericVector infectivity_curve;
  Rcpp::NumericVector susceptibility_curve;
  Rcpp::IntegerVector possible_ages;
  Rcpp::NumericVector age_distribution;

public:
  std::default_random_engine generator;

  // Constructors
  HostPopulation();
  HostPopulation(int initialS, int initialI, int initialR, double iniDay, double _contactRate, double _mu, 
                 double _wane, int seed_variant);
  HostPopulation(int initialS, int initialI, int initialR, double iniDay, double _contactRate, double _mu, 
                 double _wane, int seed_variant,
                 Rcpp::NumericVector _age_distribution, 
                 Rcpp::NumericVector _infectivity_curve, Rcpp::NumericVector _susceptibility_curve);
  HostPopulation(int initialS, int initialI, int initialR, double iniDay, double _contactRate, double _mu, double _wane, 
                 double _max_vacc_immunity, double _inf_wane_rate, double _vacc_wane_rate, int seed_variant);

  // Destructor
  ~HostPopulation();

  // Manage population temporal dynamics
  double stepForward(double new_day);
  void grow();
  void vaccinations();
  void decline();
  void seed(int variant, int number); // Function to seed new infections of a certain variant
  double contact();
  void recoveries();
  void onsets();
  void waning();
  void updateCompartments();
  
  
  // Get properties of HostPopulation
  double getContactRate();
  double getDay();
  int countSusceptibles();
  int countInfecteds();
  int countRecovereds();
  int countN();
  
  // Change parameters
  void set_contactRate(double _contactRate);
  
  // Access parameters
  double getMaxVaccImmunity();
  double getInfWaneRate();
  double getVaccWaneRate();
  
  double getInfectivity(int age);
  double getSusceptibility(int age);

  // Print out current population status
  bool inVirusVector(Virus* V);
  
  void printStatus();
  void writeHosts(std::ofstream& output, std::string filename);
  void writeViruses(std::ofstream& output, std::string filename);
  
  void testingRound(std::ofstream& output, int n_samples);

};

#endif
