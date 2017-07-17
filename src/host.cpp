#include "host.hpp"

using namespace std;

int Host::_meanBoost = 10;
int Host::_maxTitre = 10;

void Host::changeMeanBoost(int newBoost){
  _meanBoost = newBoost;
}

void Host::set_maxTitre(int newTitre){
  _maxTitre = newTitre;
}

Host::Host(){
  state = Susceptible;
  currentInfection = NULL;
  hostK = 0;
}

Host::Host(State _state, HostPopulation* _popn){
  state = _state;
  currentInfection = NULL;
  popn = _popn;
  hostK = 0;
}


Host::Host(State _state, HostPopulation* _popn, int _k){
  state = _state;
  currentInfection = NULL;
  popn = _popn;
  hostK = (double)_k;
}

Host::Host(State _state, HostPopulation* _popn, int _k, Virus* _firstInf){
state = _state;
  currentInfection = NULL;
  popn = _popn;
  hostK = (double)_k;
  infectionHistory.push_back(_firstInf);
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

int Host::decaying_boost(){
  int boost;
  double scale = (-(double)_meanBoost/(double)_maxTitre)*(double)hostK + (double)_meanBoost;
  boost = R::rpois(scale);
  if(boost < 0) boost = 0;
  return(boost);
}

void Host::recover(int cur_t){
  state = Recovered;
  int boost = decaying_boost();
  hostK += boost;

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

double Host::calculateBeta(){
  return(popn->getContactRate()*currentInfection->calculateRho(this));
}

double Host::get_hostK(){
  return(hostK);
}

default_random_engine Host::get_generator(){
  return(popn->generator);
}
