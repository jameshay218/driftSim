#include "host.hpp"

using namespace std;

Host::Host(){
  state = Susceptible;
  currentInfection=NULL;
  symptomatic=false;
}

Host::Host(State _state, HostPopulation* _popn){
  state = _state;
  symptomatic=false;
  currentInfection=NULL;
  popn = _popn;
}


Host::Host(State _state, HostPopulation* _popn, bool _symptomatic){
  symptomatic = _symptomatic;
  currentInfection=NULL;
  state = _state;
  popn = _popn;
}

Host::Host(State _state, HostPopulation* _popn, bool _symptomatic, Virus* firstInf){
  symptomatic = _symptomatic;
  currentInfection=NULL;
  state = _state;
  popn = _popn;
  infectionHistory.push_back(firstInf);
}

Host::~Host(){
  int j = infectionHistory.size();
  for(int i = 0; i < j;++i){
    if(infectionHistory[i] != NULL && infectionHistory[i] != popn->getSeedVirus()){
      delete infectionHistory[i];
    }
  }
}

bool Host::isInfected(){
  return(state == 1);
}

bool Host::isDead(){
  return(state==3);
}

bool Host::isRecovered(){
  return(state==2);
}

bool Host::isSusceptible(){
  return(state==0);
}

void Host::infect(Virus* newInfection, int cur_t){
  state = Infected;
  
  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = newInfection;
}


void Host::recover(int cur_t){
  state = Recovered;
  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = NULL;
}

Virus* Host::getCurrentVirus(){
  return(currentInfection);
}

std::vector<Virus*> Host::getInfectionHistory(){
  return(infectionHistory);
}

State Host::getState(){
  return(state);
}

void Host::addInfection(Virus* infection){
  infectionHistory.push_back(infection);
}

void Host::die(int cur_t){
  state = Dead;
  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = NULL;
}

void Host::wane(){
  state = Susceptible;
}

double Host::calculateInfectiousness(int cur_t){
  double vl = 0;
  if(currentInfection != NULL){
      vl = currentInfection->calculateViralLoad(cur_t);
  }
  return vl;
}
double Host::calculateSusceptibility(Virus* infecting_virus){
  
}

default_random_engine Host::get_generator(){
  return(popn->generator);
}
