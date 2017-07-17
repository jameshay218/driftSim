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
  std::vector<Virus*> infectionHistory;
  Virus* currentInfection;
  int hostK;

public:
  static int _meanBoost;
  static int _maxTitre;

  HostPopulation* popn;
  // Constructors
  Host();
  Host(State _state, HostPopulation* _popn, int _k, Virus* firstInf);
  Host(State _state, HostPopulation* _popn);
  Host(State _state, HostPopulation* _popn, int _k);
  ~Host();

  // Calculations/Events
  double calculateBeta();
  void infect(Virus* newInfection, int cur_t);
  void recover(int cur_t);
  void addInfection(Virus* infection);
  void wane();
  void die(int cur_t);
  double get_hostK();
  int decaying_boost();

  bool isInfected();
  bool isDead();
  bool isSusceptible();
  bool isRecovered();

  // Attribute Access
  Virus* getCurrentVirus();
  std::vector<Virus*> getInfectionHistory();
  State getState();
  std::default_random_engine get_generator();

  static void changeMeanBoost(int newBoost);
  static void set_maxTitre(int newTitre);

};

#endif
