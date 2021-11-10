#include "hostpopulation.hpp"
#include "host.hpp"
#include "virus.hpp"

using namespace std;

HostPopulation::HostPopulation(){
  day = 1;
  contactRate = 0.5;
  mu = 1.0/(70.0*365.0);
  for(int i = 0; i < 70000; ++i){
    Host* firstH = new Host();
    susceptibles.push_back(firstH);
  }
}

HostPopulation::HostPopulation(int initialS, int initialI, int initialR, double iniDay, 
                               double _contactRate, double _mu, double _wane,
                               int seed_variant){
  
  // Pointers for creating hosts and viruses
  Host* H;
  Virus* V;
  Virus* parentV;
  
  day = iniDay;
  contactRate = _contactRate;
  mu = _mu;
  wane = _wane;
  
  vacc_wane_rate=-0.03;
  inf_wane_rate=0;
  max_vacc_immunity=0.7;
  vacc_prob=0.01;
  
  Rcpp::Rcout << "Creating host population" << std::endl;
  
  // Create seed virus
  parentV = new Virus(-1, seed_variant, NULL, NULL, 0);
  seed_viruses.push_back(parentV);
  
  for(int i = 0; i < initialS;++i){
    H = new Host(Susceptible, this);
    susceptibles.push_back(H);
  }
  for(int i =0; i < initialI; ++i){
    H = new Host(Infected, this);
      // Create new virus with parameters
    V = new Virus(day, seed_variant, parentV, H, 0);
    H->infect(V,day);
    infecteds.push_back(H);
  }
  for(int i = 0; i < initialR; ++i){
    H = new Host(Recovered, this);
    recovereds.push_back(H);
  }
}

HostPopulation::HostPopulation(int initialS, int initialI, int initialR, double iniDay, 
                               double _contactRate, double _mu, double _wane,
                               double _max_vacc_immunity, double _inf_wane_rate, double _vacc_wane_rate,
                               int seed_variant){
  
  // Pointers for creating hosts and viruses
  Host* H;
  Virus* V;
  Virus* parentV;
  
  day = iniDay;
  contactRate = _contactRate;
  mu = _mu;
  wane = _wane;
  
  vacc_wane_rate=_vacc_wane_rate;
  inf_wane_rate=_inf_wane_rate;
  max_vacc_immunity=_max_vacc_immunity;
  vacc_prob=0.1;
  
  Rcpp::Rcout << "Creating host population" << std::endl;
  
  // Create seed virus
  parentV = new Virus(-1, seed_variant, NULL, NULL, 0);
  seed_viruses.push_back(parentV);
  
  for(int i = 0; i < initialS;++i){
    H = new Host(Susceptible, this);
    susceptibles.push_back(H);
  }
  for(int i =0; i < initialI; ++i){
    H = new Host(Infected, this);
    // Create new virus with parameters
    V = new Virus(day, seed_variant, parentV, H, 0);
    H->infect(V,day);
    infecteds.push_back(H);
  }
  for(int i = 0; i < initialR; ++i){
    H = new Host(Recovered, this);
    recovereds.push_back(H);
  }
}

HostPopulation::~HostPopulation(){
  int j = susceptibles.size();
  for(int i = 0; i < j; ++i){
    delete susceptibles[i];
  }
  j = infecteds.size();
  for(int i = 0; i < j; ++i){
    delete infecteds[i];
  }
  j = recovereds.size();
  for(int i = 0; i < j; ++i){
    delete recovereds[i];
  }
  j = dead.size();
  for(int i = 0; i < j; ++i){
    delete dead[i];
  }
  // Delete seed viruses
  j = seed_viruses.size();
  for(int i = 0; i < j; ++i){
      if(seed_viruses[i] != NULL){
        delete seed_viruses[i];
      }
  }
}

// Returns the number of new infections so we can calculate incidence
double HostPopulation::stepForward(double new_day){
  double new_infected;
  
  // Current day is changed
  day = new_day;

  //New births of susceptibles
  // Rcpp::Rcout << "Grow" << endl;
  grow();
  
  // Vaccinations
  vaccinations();

  // Infected population grows from transmission
  //Rcpp::Rcout << "Contact" << endl;
  new_infected = contact();

  // Infected population declines from recovery
  // Rcpp::Rcout << "Recovered" << endl;
  recoveries();

  onsets();
  
  // Some immunity wanes
  // Rcpp::Rcout << "Waning" << endl;
  waning();

  // Deaths from all compartments - MUST BE LAST EVENT
  //Rcpp::Rcout << "Decline" << endl;
  decline();
  
  //Rcpp::Rcout << "Update" << endl;
  updateCompartments();
  
  return(new_infected);
}

void HostPopulation::seed(int variant, int number){
    
    // Create new seed virus
    Virus* newSeedV = new Virus(-1, variant, NULL, NULL, 0);
    seed_viruses.push_back(newSeedV);
    
    int index;
    
    // For each seed infection, create new infection from this seed virus
    for(int i = 0; i < number; ++i){
        // Choose a random susceptible -- NOTE, WE MIGHT ACCIDENTALLY CHOOSE THE SAME PERSON TWICE
        index = floor(R::unif_rand()*(countSusceptibles()));
        
        // If this person is indeed susceptible, infect when with a child of the new variant
        if(susceptibles[index]->isSusceptible()){
            Virus* newV = new Virus(day, variant, newSeedV, susceptibles[index], susceptibles[index]->getInfectionHistory().size());
            susceptibles[index]->infect(newV,day);
            new_infecteds.push_back(susceptibles[index]); // Record new infected to add later
        }
    }
    
    // Update these new infecteds and lost susceptibles
    updateCompartments();
}

void HostPopulation::grow(){
  /*  poisson_distribution<int> poisson(mu*countN());
      int newBirths = poisson(generator);*/
  int newBirths = R::rpois(mu*countN());
  
  for(int i = 0; i < newBirths; ++i){
    Host* h = new Host(Susceptible, this);
    new_births.push_back(h);
  }
}

void HostPopulation::vaccinations(){
  // Put all hosts into one vector
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds};
  
  // Go through all hosts
  int end = tmp.size();
  int tot = 0;
  for(int i = 0; i < end; ++i){
    // Go through each State type
    tot = tmp[i].size();
    for(int j = 0; j < tot;++j){
      if(!tmp[i][j]->isVaccinated()){
        if(R::unif_rand() < vacc_prob){
          tmp[i][j]->vaccinate(day);
        }
      }
    }
  }
}


void HostPopulation::decline(){
  int newDeaths = R::rpois(mu*countSusceptibles());
  
  int totDeaths = 0;
  int index;
  for(int i = 0; i < newDeaths; ++i){
    if(countSusceptibles() > 0){
      
      index = floor(R::unif_rand()*(countSusceptibles()));
      
      dead.push_back(susceptibles[index]);
      susceptibles[index]->die(day);
      susceptibles[index] = susceptibles.back();
      susceptibles.pop_back();
      totDeaths++;
    }
  }
  
  newDeaths = R::rpois(mu*countInfecteds());
  
  for(int i = 0; i < newDeaths; ++i){
    if(countInfecteds() > 0){
      
      index = floor(R::unif_rand()*(countInfecteds()));
      
      dead.push_back(infecteds[index]);
      infecteds[index]->die(day);
      infecteds[index] = infecteds.back();
      infecteds.pop_back();
      totDeaths++;
    }
  }
  
  newDeaths = R::rpois(mu*countRecovereds());
  
  for(int i = 0; i < newDeaths; ++i){
    if(countRecovereds() > 0){
 
      index = floor(R::unif_rand()*(countRecovereds()));
  
      dead.push_back(recovereds[index]);
      recovereds[index]->die(day);
      recovereds[index] = recovereds.back();
      recovereds.pop_back();
      totDeaths++;
    }
  }
}

double HostPopulation::contact(){
  // Generate number of contacts between infecteds and susceptibles
  int totalContacts = R::rpois(contactRate*countInfecteds()*countSusceptibles()/countN());
  
  int index1 = 0;
  int index2 = 0;
  double tmp = 0;
  double number_success = 0;
  
  double probTransmission;
  double probInfection;
  
  /* For each contact, get a random I and random S. With a probability proportional to virus survival, infect the susceptible. 
     Note that the susceptible may contact multiple infecteds */
  for(int i = 0; i < totalContacts; ++i){
    if(countInfecteds() > 0){
      
      index1 = floor(R::unif_rand()*(countInfecteds()));
      index2 = floor(R::unif_rand()*(countSusceptibles()));
      tmp = R::unif_rand();
      
      // Check if given virus can infect given host
      if(susceptibles[index2]->isSusceptible()){
          // Calculate probability of transmission
          probTransmission = infecteds[index1]->calculateInfectiousness(day); // Infector has certain infectiousness
          probInfection = susceptibles[index2]->calculateSusceptibility(infecteds[index1]->getCurrentVirus(), day); // Infectee has immunity against this virus
            // Successful if transmission occurs and this gets through immunity  
          if(tmp <= probTransmission*probInfection){
              // If successful infection, use the infecting virus as the parent. 
              Virus* newV = new Virus(day, infecteds[index1]->getCurrentVirus()->getVariant(), 
                                          infecteds[index1]->getCurrentVirus(), susceptibles[index2],
                                          susceptibles[index2]->getInfectionHistory().size());
              susceptibles[index2]->infect(newV,day);
              new_infecteds.push_back(susceptibles[index2]); // Record new infected to add later
              number_success++;
          }
      }
    }
  }
  return(number_success);
}

void HostPopulation::recoveries(){
  // Go through all currently infected and check if they've recovered yet
  int n_infected = countInfecteds();
  for(int i = 0; i < n_infected; ++i){
      // Check if individual has recovered today
      if(infecteds[i]->recover(day)){
        // If so, add to recovered list and remove from infected
        new_recovereds.push_back(infecteds[i]);
        infecteds[i] = infecteds.back();
        infecteds.pop_back();
      }
  }
  // Add the newly recovered individuals back to list of infecteds so that other events can be calculated correctly.
  infecteds.insert(infecteds.end(),new_recovereds.begin(),new_recovereds.end());
}

void HostPopulation::onsets(){
    // Go through all currently infected and check if they've had symptom onset today
    int n_infected = countInfecteds();
    for(int i = 0; i < n_infected; ++i){
        // Check if individual has recovered today
        infecteds[i]->onset(day);
    }
  
  // Add the newly recovered individuals back to list of infecteds so that other events can be calculated correctly.
  infecteds.insert(infecteds.end(),new_recovereds.begin(),new_recovereds.end());
}

void HostPopulation::waning(){
    // Some individuals in recovered compartment become susceptible again
  int noWane = R::rpois(wane*countRecovereds());
  
  int index = 0;
  for(int i = 0; i < noWane; ++i){
    if(countRecovereds() > 0){
      
      index = floor(R::unif_rand()*(countRecovereds()));
      
      new_susceptibles.push_back(recovereds[index]);
      recovereds[index]->wane();
      recovereds[index] = recovereds.back();
      recovereds.pop_back();
    }
  }
  // Add the newly susceptible individuals back to list of recovereds so that other events can be calculated correctly.
  recovereds.insert(recovereds.end(),new_susceptibles.begin(),new_susceptibles.end());
}

void HostPopulation::updateCompartments(){
  infecteds.insert(infecteds.end(),new_infecteds.begin(),new_infecteds.end());
  recovereds.insert(recovereds.end(),new_recovereds.begin(),new_recovereds.end());
  susceptibles.insert(susceptibles.end(),new_susceptibles.begin(),new_susceptibles.end());
  susceptibles.insert(susceptibles.end(),new_births.begin(),new_births.end());

  susceptibles.erase(remove_if(susceptibles.begin(),susceptibles.end(), [](Host *h){return(!h->isSusceptible());}),susceptibles.end());
  infecteds.erase(remove_if(infecteds.begin(),infecteds.end(), [](Host* o){return(!o->isInfected());}),infecteds.end());
  recovereds.erase(remove_if(recovereds.begin(),recovereds.end(), [](Host* o){return(!o->isRecovered());}),recovereds.end());

  new_infecteds.clear();
  new_recovereds.clear();
  new_susceptibles.clear();
  new_births.clear();
}


int HostPopulation::countSusceptibles(){
  return(susceptibles.size());
}

int HostPopulation::countInfecteds(){
  return(infecteds.size());
}

int HostPopulation::countRecovereds(){
  return(recovereds.size());
}

double HostPopulation::getDay(){
  return(day);
}

int HostPopulation::countN(){
  return(infecteds.size() + recovereds.size() + susceptibles.size());
}

double HostPopulation::getContactRate(){
  return(contactRate);
}

void HostPopulation::set_contactRate(double _contactRate){
  contactRate = _contactRate;
}


double HostPopulation::getMaxVaccImmunity(){
  return(max_vacc_immunity);
}
double HostPopulation::getInfWaneRate(){
  return(inf_wane_rate);
}
double HostPopulation::getVaccWaneRate(){
  return(vacc_wane_rate);
}


bool HostPopulation::inVirusVector(Virus* V){
  int j = seed_viruses.size();
  bool result = false;
  for(int i = 0; i < j; i++){
    if(V == seed_viruses[i]){
      result = true;
    }
  }
  return(result);  
}


void HostPopulation::printStatus(){
  Rcpp::Rcout << "Current day: " << day << endl;
  Rcpp::Rcout << "Susceptible: " << susceptibles.size() << endl;
  Rcpp::Rcout << "Infecteds: " << infecteds.size() << endl;
  Rcpp::Rcout << "Recovereds: " << recovereds.size() << endl;
  Rcpp::Rcout << "Died: " << dead.size() << endl;
  Rcpp::Rcout << "Total: " << countN() << endl;
}


/* ------------------------------------------- FILE MANIPULATION CODE ---------------------------------------------- */
 
void HostPopulation::testingRound(std::ofstream& output, int n_sample){
  // Put all hosts into one vector
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds};
  
  int N_total = countN();
  
  double prob_sample = (double) n_sample/ (double) N_total;
  
  // Go through all hosts
  int end = tmp.size();
  int tot = 0;
  for(int i = 0; i < end; ++i){
    // Go through each State type
    tot = tmp[i].size();
    for(int j = 0; j < tot;++j){
      // If sampled
      if(R::unif_rand() < prob_sample){
        // Write day
        output << day << ",";
        // Write state
        output << tmp[i][j]->getState() << ",";
        // If currently infected, write current viral load
        if(tmp[i][j]->isInfected() == 1){
          output << tmp[i][j]->getCurrentVirus()->calculateViralLoad(day) << ",";
          output << tmp[i][j]->getCurrentVirus()->getVariant() << ",";
        }
        else {
          output << "0" << ",";
          output << "NA" << ",";
        }
        output << tmp[i][j]->isSymptomatic() << endl;
      }
    }
  }
}

void HostPopulation::writeHosts(std::ofstream& output, std::string filename){
  Rcpp::Rcout << "#########################" << endl;
  Rcpp::Rcout << "Writing hosts to csv..." << endl;
  Rcpp::Rcout << "Location: " << filename << endl;
  output.open(filename);
  
  // Put all hosts into one vector
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds};

  output << "state,last_vid,cur_inf,symptomatic,vaccine_time,inf_history" << endl;
  
  // Go through all hosts
  int end = tmp.size();
  int tot = 0;
  int inf_hist_size=0;
  for(int i = 0; i < end; ++i){
    // Go through each State type
    tot = tmp[i].size();
    for(int j = 0; j < tot;++j){
        // Write current state
        output << tmp[i][j]->getState() << ",";
            
        // Write ID of latest infection
        if(tmp[i][j]->getInfectionHistory().size() > 0){
            output << tmp[i][j]->getInfectionHistory()[0]->getId() << ",";
        }
        else {
            output << 0 << ",";
        }
        
        // If currently infected, write current infection ID
        if(tmp[i][j]->isInfected() == 1){
    	    output << tmp[i][j]->getCurrentVirus()->getId() << ",";
        }
        else {
            output << tmp[i][j]->isInfected() << ",";
        }
        output << tmp[i][j]->isSymptomatic() << ",";
        output << tmp[i][j]->getVaccTime() << ",";
        inf_hist_size =tmp[i][j]->getInfectionHistory().size();
        if(inf_hist_size > 0){
          for(int k = 0; k < inf_hist_size; ++k){
            output << tmp[i][j]->getInfectionHistory()[k]->getId() << "-";
          }
        }
        output << std::endl;
    }
  }
  output.close();
  Rcpp::Rcout << "Hosts writing complete" << endl;
  Rcpp::Rcout << "#########################" << endl << endl;
}

void HostPopulation::writeViruses(std::ofstream& output, std::string filename){
  Rcpp::Rcout << "#########################" << endl;
  Rcpp::Rcout << "Writing viruses to csv..." << endl;
  output.open(filename);
  int x = 0;
  
  // Get all hosts from the simulation
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds, dead};
  vector<Virus*> viruses;
  int q = tmp.size();  
  int j = 0;
  
  // Go through each Host state
  for(int y = 0; y < q; ++y){
        j = tmp[y].size();
        // Within this state, get all hosts one at a time
        for(int i = 0; i < j; ++i){
            // If host is currently infected, add this virus to list to be printed
          if(tmp[y][i]->isInfected()){
            if(!inVirusVector(tmp[y][i]->getCurrentVirus())){
    	        viruses.push_back(tmp[y][i]->getCurrentVirus());
            }
          }
          // Then, go through the infection history and add each other virus
          x = tmp[y][i]->getInfectionHistory().size();
          if(x > 0){
    	    for(int ii = 0; ii < x; ++ii){
    	      if(!inVirusVector(tmp[y][i]->getInfectionHistory()[ii])){
    	        viruses.push_back(tmp[y][i]->getInfectionHistory()[ii]);
    	      }
    	    }
          }
        }
  }
  // Sort viruses into order before printing
  Rcpp::Rcout << "Sorting viruses..." << endl;
  sort(viruses.begin(),viruses.end(),[](Virus* lhs, Virus* rhs){
      return(lhs->getId() < rhs->getId());
    });
  Rcpp::Rcout << "Writing output..." << endl;
  
  output << "vid,variant_id,variant,birth,death,parentid,parent_age_at_creation,infectionNo,tg,tp,to,tw,alpha,symptomatic,infectiousnessMax,infectiousnessGradient,infectiousnessInflection,host_vacc_status" << endl;
  j = seed_viruses.size();
  
  Rcpp::Rcout << "...writing seed viruses..." << endl;
  Rcpp::Rcout << "Writing viruses complete" << endl;
  Rcpp::Rcout << "#########################" << endl << endl;

  if(j > 0){
      for(int i = 0; i < j; i++){
        output << seed_viruses[i]->getId() << ",";
        output << seed_viruses[i]->getVariantId() << ",";
        output << seed_viruses[i]->getVariant() << ",";
        output << seed_viruses[i]->getBirth() << ",";
        output << seed_viruses[i]->getDeath() << ",";
        if(seed_viruses[i]->getParent() != NULL){
	        output << seed_viruses[i]->getParent()->getId() << ",";
          output << seed_viruses[i]->getAgeOfParentAtBirth() << ",";
        } else {
	        output << "0" << ",";
          output << "0" << ",";
        }
        output << seed_viruses[i]->getInfectionNo() << ",";
        
        // Kinetics pars
        output << seed_viruses[i]->get_tg() << ",";
        output << seed_viruses[i]->get_tp() << ",";
        output << seed_viruses[i]->get_to() << ",";
        output << seed_viruses[i]->get_tw() << ",";
        output << seed_viruses[i]->get_alpha() << ",";
        output << seed_viruses[i]->get_symptomatic() << ",";
        output << seed_viruses[i]->get_infectiousnessMax() << ",";
        output << seed_viruses[i]->get_infectiousnessGradient() << ",";
        output << seed_viruses[i]->get_infectiousnessInflection() << ",";
        output << seed_viruses[i]->get_host_vacc();
        
        
        output << endl;
      }
  }
  Rcpp::Rcout << "Seed viruses saved" << endl;
  Rcpp::Rcout << "...writing all viruses..." << endl;

  j = viruses.size();
  for(int i = 0; i < j;++i){
      output << viruses[i]->getId() << ",";
      output << viruses[i]->getVariantId() << ",";
      output << viruses[i]->getVariant() << ",";
      output << viruses[i]->getBirth() << ",";
      output << viruses[i]->getDeath() << ",";
      if(viruses[i]->getParent() != NULL){
          output << viruses[i]->getParent()->getId() << ",";
      } else {
          output << "0" << ",";
      }
      output << viruses[i]->getAgeOfParentAtBirth() << ",";
      output << viruses[i]->getInfectionNo() << ",";
      
      // Kinetics pars
      output << viruses[i]->get_tg() << ",";
      output << viruses[i]->get_tp() << ",";
      output << viruses[i]->get_to() << ",";
      output << viruses[i]->get_tw() << ",";
      output << viruses[i]->get_alpha() << ",";
      output << viruses[i]->get_symptomatic() << ",";
      output << viruses[i]->get_infectiousnessMax() << ",";
      output << viruses[i]->get_infectiousnessGradient() << ",";
      output << viruses[i]->get_infectiousnessInflection() << ",";
      output << viruses[i]->get_host_vacc();
      
      output << endl;
  }
  output.close();
  Rcpp::Rcout << "Writing viruses complete" << endl;
  Rcpp::Rcout << "#########################" << endl << endl;
 }
 

