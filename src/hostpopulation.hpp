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

public:
  std::default_random_engine generator;

  // Constructors
  HostPopulation();
  HostPopulation(int initialS, int initialI, int initialR, double iniDay, double _contactRate, double _mu, double _wane, int seed_variant);

  // Destructor
  ~HostPopulation();

  // Manage population temporal dynamics
  double stepForward(double new_day);
  void grow();
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
  
  // CHange parameters
  void set_contactRate(double _contactRate);

  // Print out current population status
  bool inVirusVector(Virus* V);
  
  void printStatus();
  void writeHosts(std::ofstream& output, std::string filename);
  void writeViruses(std::ofstream& output, std::string filename);
  
  void testingRound(std::ofstream& output, int n_samples);

};

#endif
