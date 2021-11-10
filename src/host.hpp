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
  
  // Vaccination and symptom status
  bool symptomatic;
  bool vaccinated; 
  double vaccination_time; // Time of vaccination in days
  
  int age; // Age in years

public:

  HostPopulation* popn;
  
  // Constructors
  Host();
  Host(State _state, HostPopulation* _popn); // Constructor for susceptible
  Host(State _state, HostPopulation* _popn, int _age); // Constructor for susceptible with age
  Host(State _state, HostPopulation* _popn, Virus* firstInf); // Constructor if currently infected or recovered with one past infection
  Host(State _state, HostPopulation* _popn, Virus* firstInf, int _age); // Constructor if currently infected or recovered with one past infection + age
  ~Host();

  // Calculations
  double calculateInfectiousness(double cur_t); // Get infectiousness based on viral kinetics
  double calculateSusceptibility(Virus* infecting_virus, double cur_t); // Get susceptibility based on infecting virus and infection history
  
  // Events
  void infect(Virus* newInfection, double cur_t);
  void addInfection(Virus* infection);
  bool recover(double cur_t);
  void onset(double cur_t);
  void wane();
  void die(double cur_t);
  
  void vaccinate(double cur_t);

  // Attribute Access
  Virus* getCurrentVirus();
  std::vector<Virus*> getInfectionHistory();
  State getState();
  
  bool isInfected();
  bool isDead();
  bool isSusceptible();
  bool isRecovered();
  bool isSymptomatic();
  bool isVaccinated();
  double getVaccTime();
  int getAge();
  
};

#endif
