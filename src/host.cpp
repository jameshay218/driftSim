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


Host::Host(State _state, HostPopulation* _popn, Virus* firstInf){
  symptomatic = false;
  currentInfection=NULL;
  state = _state;
  popn = _popn;
  infectionHistory.push_back(firstInf);
}

Host::~Host(){
  int j = infectionHistory.size();
  for(int i = 0; i < j;++i){
    if(infectionHistory[i] != NULL){
      delete infectionHistory[i];
    }
  }
}

bool Host::isSusceptible(){
  return(state==0);
}

bool Host::isInfected(){
  return(state == 1);
}

bool Host::isRecovered(){
  return(state==2);
}

bool Host::isDead(){
  return(state==3);
}
bool Host::isSymptomatic(){
  return(symptomatic);
}

void Host::infect(Virus* newInfection, int cur_t){
  state = Infected;
  
  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = newInfection;
}

// If has current infection and viral load has waned to 0, then recover
bool Host::recover(int cur_t){
  // Only recover if currently infected and enough time has elapsed to recover
  if(currentInfection != NULL && currentInfection->hasRecovered(cur_t)){
    state = Recovered;
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
    currentInfection = NULL;
    symptomatic=false;
    return true;
  }
  return false;
}

// Check if symptom onset is today. If so, set state to symptomatic
void Host::onset(int cur_t){
  if(currentInfection != NULL && currentInfection->hasOnset(cur_t)){
    symptomatic=true;
  }
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

// If have a current infection, then get infectiousness of this virus
double Host::calculateInfectiousness(int cur_t){
  double infectiousness = 0;
  if(currentInfection != NULL){
    infectiousness = currentInfection->getInfectiousness(cur_t);
  }
  return infectiousness;
}


// Given an attempted infecting virus, find maximum cross immunity in infection history
// Can replace this with something like titer-based immunity in the future
double Host::calculateSusceptibility(Virus* infectingVirus){
  int infhist_size = infectionHistory.size();
  
  // Assume for now no immunity
  double max_immunity = 0;
  double tmp_immunity = 0;
  
  // If have past infections, loop through and find maximum cross immunity with infecting virus
  if(infhist_size > 0) {
    for(int i = 0; i < infhist_size; ++i){
      tmp_immunity = Virus::getCrossImmunity(infectingVirus, infectionHistory[i]);
      if(tmp_immunity > max_immunity){
        max_immunity = tmp_immunity;
      }
    }
  }
  return(1.0 - max_immunity);
}
