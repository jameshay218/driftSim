#include "host.hpp"

using namespace std;

Host::Host(){
  state = Susceptible;
  currentInfection=NULL;
  symptomatic=false;
  
  vaccinated=false;
  vaccination_time=-1;
  age=0;
;}

Host::Host(State _state, HostPopulation* _popn){
  state = _state;
  symptomatic=false;
  currentInfection=NULL;
  popn = _popn;
  
  vaccinated=false;
  vaccination_time=-1;
  age=0;
}

Host::Host(State _state, HostPopulation* _popn, int _age){
  state = _state;
  symptomatic=false;
  currentInfection=NULL;
  popn = _popn;
  
  vaccinated=false;
  vaccination_time=-1;
  age=_age;
}


Host::Host(State _state, HostPopulation* _popn, Virus* firstInf){
  symptomatic = false;
  currentInfection=NULL;
  state = _state;
  popn = _popn;
  infectionHistory.push_back(firstInf);
  
  vaccinated=false;
  vaccination_time=-1;
  age=0;
}

Host::Host(State _state, HostPopulation* _popn, Virus* firstInf, int _age){
  symptomatic = false;
  currentInfection=NULL;
  state = _state;
  popn = _popn;
  infectionHistory.push_back(firstInf);
  
  vaccinated=false;
  vaccination_time=-1;
  age=_age;
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
  return(state==Susceptible);
}

bool Host::isInfected(){
  return(state==Infected);
}

bool Host::isRecovered(){
  return(state==Recovered);
}

bool Host::isDead(){
  return(state==Dead);
}
bool Host::isSymptomatic(){
  return(symptomatic);
}
bool Host::isVaccinated(){
  return vaccinated;
}

void Host::infect(Virus* newInfection, double cur_t){
  state = Infected;
  
  if(currentInfection != NULL){
    currentInfection->kill(cur_t);
    infectionHistory.push_back(currentInfection);
  }
  currentInfection = newInfection;
}

void Host::vaccinate(double cur_t){
  vaccinated=true;
  vaccination_time=cur_t;
}

// If has current infection and viral load has waned to 0, then recover
bool Host::recover(double cur_t){
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
void Host::onset(double cur_t){
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

double Host::getVaccTime(){
  return(vaccination_time);
}

int Host::getAge(){
  return(age);
}

void Host::addInfection(Virus* infection){
  infectionHistory.push_back(infection);
}

void Host::die(double cur_t){
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
double Host::calculateInfectiousness(double cur_t){
  double infectiousness = 0;
  if(currentInfection != NULL){
    //Rcpp::Rcout << "Age: " << age << "; infectiousness: " << popn->getInfectivity(age) << std::endl;
    infectiousness = currentInfection->getInfectiousness(cur_t)*popn->getInfectivity(age);
  }
  return infectiousness;
}


// Given an attempted infecting virus, find maximum cross immunity in infection history
// Can replace this with something like titer-based immunity in the future
double Host::calculateSusceptibility(Virus* infectingVirus, double cur_t){
  int infhist_size = infectionHistory.size();
  
  // Assume for now no immunity
  double max_immunity = 0;
  double tmp_immunity = 0;
  double infection_time = 0;
  
  // If have past infections, loop through and find maximum cross immunity with infecting virus
  if(infhist_size > 0) {
    for(int i = 0; i < infhist_size; ++i){
      infection_time = infectionHistory[i]->getDeath();
      tmp_immunity = Virus::getCrossImmunity(infectingVirus, infectionHistory[i]);
      tmp_immunity = tmp_immunity*exp(popn->getInfWaneRate()*(cur_t-infection_time));
      if(tmp_immunity > max_immunity){
        max_immunity = tmp_immunity;
      }
    }
  }
  
  if(vaccinated){
    tmp_immunity = popn->getMaxVaccImmunity()*exp(popn->getVaccWaneRate()*(cur_t-vaccination_time));
    if(tmp_immunity > max_immunity){
      max_immunity = tmp_immunity;
    }
  }
  
  return((1.0 - max_immunity)*popn->getSusceptibility(age));
}
