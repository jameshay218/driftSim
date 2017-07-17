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

  int day;
  double contactRate;
  double mu;
  double wane;
  double gamma;

  Virus* seedVirus;

public:
  std::default_random_engine generator;

  // Constructors
  HostPopulation();
  HostPopulation(int initialS, int initialI, int initialR, int iniDay, double _contactRate, double _mu, double _wane, double _gamma, double _iniBindingAvid, double initialDistance);
  HostPopulation(int initialS, int initialI, int initialR, int iniDay, double _contactRate, double _mu, double _wane, double _gamma, double _iniBindingAvid, double initialDistance, Rcpp::NumericVector startingK);
  ~HostPopulation();

  // Manage population temporal dynamics
  void stepForward(int new_day);
  void grow();
  void decline();
  void contact();
  void recoveries();
  void waning();
  void mutations();
  void updateCompartments();

  // Get properties of HostPopulation
  double getContactRate();
  int getDay();
  int countSusceptibles();
  int countInfecteds();
  int countRecovereds();
  int countN();

  Virus* getSeedVirus();

  // Print out current population status
  void printStatus();
  void writeHosts(std::ofstream& output, std::string filename);
  void readHosts(std::string hostFilename, std::string virusFilename);
  std::map<int, Virus*> readViruses(std::string virusFilename);
  void writeViruses(std::ofstream& output, std::string filename, bool savingState);
  void virusPairwiseMatrix(std::ofstream& output, std::string filename, int sampSize);
  
  Rcpp::NumericVector getHostKDist();

};

#endif
