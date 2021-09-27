#ifndef HOST_HPP
#define HOST_HPP

#include <vector>
#include<cstdlib>

#include "hostpopulation.hpp"
#include "virus.hpp"

enum State {Susceptible, Infected, Recovered, Dead};

class Host{
private:
  State state;
  
  std::vector<Virus*> infectionHistory; // Keep a vector of past infections
  Virus* currentInfection; // Pointer to currently infecting virus
  
  bool symptomatic; // Symptomatic or not?
  
  int infection_time; 
  

public:

  HostPopulation* popn;
  
  // Constructors
  Host();
  Host(State _state, HostPopulation* _popn); // Constructor for susceptible
  Host(State _state, HostPopulation* _popn, Virus* firstInf); // Constructor if currently infected or recovered with one past infection
  ~Host();

  // Calculations/Events
  double calculateInfectiousness(int cur_t); // Get infectiousness based on viral kinetics
  double calculateSusceptibility(Virus* infecting_virus); // Get susceptibility based on infecting virus and infection history
  void infect(Virus* newInfection, int cur_t);
  void addInfection(Virus* infection);
  void recover(int cur_t);
  void wane();
  void die(int cur_t);

  bool isInfected();
  bool isDead();
  bool isSusceptible();
  bool isRecovered();

  // Attribute Access
  Virus* getCurrentVirus();
  std::vector<Virus*> getInfectionHistory();
  State getState();
  std::default_random_engine get_generator();
  
};

#endif
