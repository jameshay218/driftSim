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

HostPopulation::HostPopulation(int initialS, int initialI, int initialR, int iniDay, 
                               double _contactRate, double _mu, double _wane, double _gamma){
  
  // Pointers for creating hosts and viruses
  Host* H;
  Virus* V;
  Virus* parentV;
  
  day = iniDay;
  contactRate = _contactRate;
  mu = _mu;
  wane = _wane;
  gamma = _gamma;
  
  // Create seed virus
  parentV = new Virus(0, 0, NULL, NULL, 0, 
                      0, 0, 0, 0, 0);
  
  for(int i = 0; i < initialS;++i){
    H = new Host(Susceptible, this);
    susceptibles.push_back(H);
  }
  for(int i =0; i < initialI; ++i){
    H = new Host(Infected, this);
      // Create new virus with parameters
    V = new Virus(day, 0, parentV, H, 0);
    H->infect(day);
    infecteds.push_back(H);
  }
  for(int i = 0; i < initialR; ++i){
    H = new Host(Recovered, this);
    recovereds.push_back(H);
  }
}


HostPopulation::HostPopulation(int initialS, int initialI, int initialR, int iniDay, 
                               double _contactRate, double _mu, double _wane, double _gamma){
  Host* H;
  int hostIndex = 0;
  day = iniDay;
  contactRate = _contactRate;
  mu = _mu;
  wane = _wane;
  gamma = _gamma;

  for(int i = 0; i < initialS;++i){
    H = new Host(Susceptible, this);
    susceptibles.push_back(H);
  }
  for(int i =0; i < initialI; ++i){
    H = new Host(Infected, this);
    H->infect(day);
    infecteds.push_back(H);
  }
  for(int i = 0; i < initialR; ++i){
    H = new Host(Recovered, this);
    H->addInfection();
    recovereds.push_back(H);
  }
}

HostPopulation::~HostPopulation(){
  // Rcpp::Rcout << "Host Population Destructor" << endl;
  int j = susceptibles.size();
  for(int i = 0; i < j; ++i){
    //Rcpp::Rcout << i << endl;
    //    Rcpp::Rcout << susceptibles[i]->getState() << endl;
    delete susceptibles[i];
  }
  //Rcpp::Rcout << "S deleted" << endl;
  j = infecteds.size();
  for(int i = 0; i < j; ++i){
    delete infecteds[i];
  }
  // Rcpp::Rcout << "I deleted" << endl;
  j = recovereds.size();
  for(int i = 0; i < j; ++i){
    delete recovereds[i];
  }
  //  Rcpp::Rcout << "R deleted" << endl;
  j = dead.size();
  for(int i = 0; i < j; ++i){
    // if(!dead[i]->isDead()){
      // Rcpp::Rcout << dead[i]->getState() << endl;
    //  }
    delete dead[i];
  }
}

// Returns the number of new infections so we can calculate incidence
double HostPopulation::stepForward(int new_day){
  double new_infected;
  
  // Current day is changed
  day = new_day;

  //New births of susceptibles
  // Rcpp::Rcout << "Grow" << endl;
  grow();

  // Infected population grows from transmission
  //Rcpp::Rcout << "Contact" << endl;
  new_infected = contact();

  // Infected population declines from recovery
  // Rcpp::Rcout << "Recovered" << endl;
  recoveries();

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


void HostPopulation::grow(){
  /*  poisson_distribution<int> poisson(mu*countN());
      int newBirths = poisson(generator);*/
  int newBirths = R::rpois(mu*countN());
  
  for(int i = 0; i < newBirths; ++i){
    Host* h = new Host(Susceptible, this);
    new_births.push_back(h);
  }
}

void HostPopulation::decline(){
  /*poisson_distribution<int> poissonS(mu*countSusceptibles());
  poisson_distribution<int> poissonI(mu*countInfecteds());
  poisson_distribution<int> poissonR(mu*countRecovereds());
  int newDeaths = poissonS(generator);*/
  
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
  
  //  newDeaths = poissonI(generator);
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
  
  //  newDeaths = poissonR(generator);
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
  double tmpDist;
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
	      // Find antigenic distance to host's immunity
	      /*
      	tmpDist = infecteds[index1]->getCurrentVirus()->findDistanceToHost(susceptibles[index2]);
      	if(tmp <= infecteds[index1]->getCurrentVirus()->calculateRho(tmpDist, susceptibles[index2])){
      	  // If successful infection, use the infecting virus as the parent. Get distance to new host's immunity.
      	  Virus* newV = new Virus(infecteds[index1]->getCurrentVirus(), susceptibles[index2], day, tmpDist, susceptibles[index2]->getInfectionHistory().size());
      	  susceptibles[index2]->infect(newV,day);
      	  new_infecteds.push_back(susceptibles[index2]); // Record new infected to add later
      	  number_success++;
      	}
	       */
      }
    }
  }
  return(number_success);
}

void HostPopulation::recoveries(){
  // Random number of recoveries of current infecteds
  /*  poisson_distribution<int> poisson(gamma*countInfecteds());
      int noRecovered = poisson(generator);*/
  
  int noRecovered = R::rpois(gamma*countInfecteds());
  

  int index = 0;
  // For each recovery, choose a random infected individual and add them to the list of newly recovered individuals
  for(int i = 0; i < noRecovered; ++i){
    if(countInfecteds() > 0){
      
      index = floor(R::unif_rand()*(countInfecteds()));
      
      // Record the new recovered individual and temporarily delete from the list of infecteds
      new_recovereds.push_back(infecteds[index]);
      infecteds[index]->recover(day);
      infecteds[index] = infecteds.back();
      infecteds.pop_back();
    }
  }
  // Add the newly recovered individuals back to list of infecteds so that other events can be calculated correctly.
  infecteds.insert(infecteds.end(),new_recovereds.begin(),new_recovereds.end());
}

void HostPopulation::waning(){
  /*  poisson_distribution<int> poisson(wane*countRecovereds());
      int noWane = poisson(generator);*/
  
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

int HostPopulation::getDay(){
  return(day);
}

int HostPopulation::countN(){
  return(infecteds.size() + recovereds.size() + susceptibles.size());
}

double HostPopulation::getContactRate(){
  return(contactRate);
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
void HostPopulation::writeHosts(std::ofstream& output, std::string filename){
  Rcpp::Rcout << "#########################" << endl;
  Rcpp::Rcout << "Writing hosts to csv..." << endl;
  Rcpp::Rcout << "Location: " << filename << endl;
  output.open(filename);
  
  // Put all hosts into one vector
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds};

  output << "state,last_vid,cur_inf,hostK" << endl;
  
  int end = tmp.size();
  int tot = 0;
  for(int i = 0; i < end; ++i){
    tot = tmp[i].size();
    for(int j = 0; j < tot;++j){
      output << tmp[i][j]->getState() << ",";
      if(tmp[i][j]->getInfectionHistory().size() > 0){
	output << tmp[i][j]->getInfectionHistory()[0]->getId() << ",";
      }
      else {
	output << 0 << ",";
      }
      if(tmp[i][j]->isInfected() == 1){
	output << tmp[i][j]->getCurrentVirus()->getId() << ",";
      }
      else {
	output << tmp[i][j]->isInfected() << ",";
      }
      output << tmp[i][j]->get_hostK() << endl;
    }
  }
  output.close();
  Rcpp::Rcout << "Hosts writing complete" << endl;
  Rcpp::Rcout << "#########################" << endl << endl;
}

void HostPopulation::writeViruses(std::ofstream& output, std::string filename, bool savingState){
  Rcpp::Rcout << "#########################" << endl;
  Rcpp::Rcout << "Writing viruses to csv..." << endl;
  output.open(filename);
  int x = 0;
  vector<vector<Host*>> tmp = {susceptibles, infecteds, recovereds, dead};
  vector<Virus*> viruses;
  int q = tmp.size();  
  int j = 0;
  for(int y = 0; y < q; ++y){
    j = tmp[y].size();
    for(int i = 0; i < j; ++i){
      if(tmp[y][i]->isInfected()){
	if(tmp[y][i]->getCurrentVirus() != seedVirus){
	  viruses.push_back(tmp[y][i]->getCurrentVirus());
	}
      }
      x = tmp[y][i]->getInfectionHistory().size();
      if(x > 0){
	for(int ii = 0; ii < x; ++ii){
	  if(tmp[y][i]->getInfectionHistory()[ii]!=seedVirus){
	    viruses.push_back(tmp[y][i]->getInfectionHistory()[ii]);
	  }
	}
      }
    }
  }
  Rcpp::Rcout << "Sorting viruses..." << endl;
  sort(viruses.begin(),viruses.end(),[](Virus* lhs, Virus* rhs){
      return(lhs->getId() < rhs->getId());
    });
  Rcpp::Rcout << "Writing output..." << endl;
  if(!savingState){
    Rcpp::Rcout << "... for phylogenetic tree" << endl;
    output << "vid, birth, death,birth_from_parent, parentid, bindingAvidityIni, bindingAvidityFinal, infectionNo, distToParent, hostImmunity, hostInfections, immJ,distHost, distRoot,changeFromV, changeFromR" << endl;
  }
  else {
    Rcpp::Rcout << "... for save state" << endl;
    output << "vid,birth,death,birth_from_parent,parentid,level,bindingAvidIni,bindingAvid,infectionNo,distToParent,distToRoot, distToHost, immJ, changeFromV, changeFromR" << endl;
  }
  if(seedVirus != NULL){
    if(!savingState){
      output << seedVirus->getId() << ",";
      output << seedVirus->getBirth() << ",";
      output << seedVirus->getDeath() << ",";
      output << seedVirus->getTimeOfInfector() << ",";
      if(seedVirus->getParent() != NULL){
	output << seedVirus->getParent()->getId() << ",";
      } else {
	output << "0" << ",";
      }
      output << seedVirus->getIniBindingAvid() << ",";
      output << seedVirus->getBindingAvid() << ",";
      output << seedVirus->getInfectionNo() << ",";
      output << seedVirus->getDistance() << ",";
      if(seedVirus->getHost() != NULL){
	output << seedVirus->getK() << ",";
	output << seedVirus->getHost()->getInfectionHistory().size() << ",";
      } else {
	output << "0" << "," << "0" << ",";
      }
      output << seedVirus->getJ() << ",";
      output << seedVirus->getDistHost() << ",";
      output << seedVirus->getDistRoot() << ",";
      output << seedVirus->getChangeFromV() << ",";
      output << seedVirus->getChangeFromR();
    }
    else {
      output << seedVirus->getId() << ",";
      output << seedVirus->getBirth() << ",";
      output << seedVirus->getDeath() << ",";
      output << seedVirus->getTimeOfInfector() << ",";
      if(seedVirus->getParent() != NULL){
	output << seedVirus->getParent()->getId() << ",";
      } else {
	output << "0" << ",";
      }
      output << seedVirus->getLevel() << ",";
      output << seedVirus->getIniBindingAvid() << ",";
      output << seedVirus->getBindingAvid() << ",";
      output << seedVirus->getInfectionNo() << ",";
      output << seedVirus->getDistance() << ",";
      output << seedVirus->getDistRoot() << ",";
      output << seedVirus->getDistHost() << ",";
      output << seedVirus->getK() << ",";
      output << seedVirus->getChangeFromV() << ",";
      output << seedVirus->getChangeFromR();
    }
    output << endl;
    Rcpp::Rcout << "Seed virus saved" << endl;
  }

  j = viruses.size();
  for(int i = 0; i < j;++i){
    if(!savingState){
      output << viruses[i]->getId() << ",";
      output << viruses[i]->getBirth() << ",";
      output << viruses[i]->getDeath() << ",";
      output << viruses[i]->getTimeOfInfector() << ",";
      if(viruses[i]->getParent() != NULL){
	output << viruses[i]->getParent()->getId() << ",";
      } else {
	output << "0" << ",";
      }
      output << viruses[i]->getIniBindingAvid() << ",";
      output << viruses[i]->getBindingAvid() << ",";
      output << viruses[i]->getInfectionNo() << ",";
      output << viruses[i]->getDistance() << ",";
      output << viruses[i]->getK() << ",";
      output << viruses[i]->getHost()->getInfectionHistory().size() << ",";
      output << viruses[i]->getJ() << ",";
      output << viruses[i]->getDistHost() << ",";
      output << viruses[i]->getDistRoot() << ",";
      output << viruses[i]->getChangeFromV() << ",";
      output << viruses[i]->getChangeFromR();
    }
    else {
      output << viruses[i]->getId() << ",";
      output << viruses[i]->getBirth() << ",";
      output << viruses[i]->getDeath() << ",";
      output << viruses[i]->getTimeOfInfector() << ",";
      if(viruses[i]->getParent() != NULL){
	output << viruses[i]->getParent()->getId() << ",";
      } else {
	output << "0" << ",";
      }
      output << viruses[i]->getLevel() << ",";
      output << viruses[i]->getIniBindingAvid() << ",";
      output << viruses[i]->getBindingAvid() << ",";
      output << viruses[i]->getInfectionNo() << ",";
      output << viruses[i]->getDistance() << ",";
      output << viruses[i]->getDistRoot() << ",";
      output << viruses[i]->getDistHost() << ",";
      output << viruses[i]->getK() << ",";
      output << viruses[i]->getChangeFromV() << ",";
      output << viruses[i]->getChangeFromR();
    }
    output << endl;
  }

  output.close();
  Rcpp::Rcout << "Writing viruses complete" << endl;
  Rcpp::Rcout << "#########################" << endl << endl;
  }

  void HostPopulation::virusPairwiseMatrix(std::ofstream& output, std::string filename, int sampSize){
    Rcpp::Rcout << "#########################" << endl;
    Rcpp::Rcout << "Writing pairwise matrix" << endl;
    vector<Host*> hosts;
    vector<Virus*> viruses;
    vector<int> sampIndices;
    output.open(filename);
    hosts.insert(hosts.end(),susceptibles.begin(),susceptibles.end());
    hosts.insert(hosts.end(),infecteds.begin(),infecteds.end());
    hosts.insert(hosts.end(),recovereds.begin(),recovereds.end());
    hosts.insert(hosts.end(),dead.begin(),dead.end());
    int j = hosts.size();
    int x = 0;
    Rcpp::Rcout << "Accumulating all viruses..." << endl;
  for(int i = 0; i < j; ++i){
    x = hosts[i]->getInfectionHistory().size();
    if(x > 0){
      for(int ii = 0; ii < x; ++ii) {
	viruses.push_back(hosts[i]->getInfectionHistory()[ii]);
      }
    }
  }
  Rcpp::Rcout << "Total number of viruses: " << viruses.size() << endl;
  j = viruses.size();
  for(int i = 0; i < j;++i){
    sampIndices.push_back(i);
  }
  random_shuffle(sampIndices.begin(),sampIndices.end());
  sampIndices.resize(sampSize);
  sort(sampIndices.begin(),sampIndices.end());  

  Rcpp::Rcout << "Sorting virus vector..." << endl;

  sort(viruses.begin(),viruses.end(),[](Virus* lhs, Virus* rhs){
      return(lhs->getId() < rhs->getId());
    });
  
  Rcpp::Rcout << "Finding pairwise matrix..." << endl;
 
  j = sampIndices.size();
  int index = 0;
  int index2 = 0;
  for(int row = 0; row < j; ++row){
    index = sampIndices[row];
    int col = 0;
    while(col < row-1){
      output << ",";
      col++;
    }
    col = row;
    while(col < j){
      output << ",";
      index2 = sampIndices[col++];
      output << Virus::getAntigenicDistance(viruses[index],viruses[index2]);
    }
    output << endl;
  }
  Rcpp::Rcout << "Pairwise output complete" << endl;
  Rcpp::Rcout << "#########################" << endl << endl;
}  
  



void HostPopulation::readHosts(std::string hostFilename, std::string virusFilename){
Rcpp::Rcout << "#########################" << endl;
 Rcpp::Rcout << "Reading hosts from csv for initial conditions..." << endl;
  int tmpSize = 0;
  ifstream hostFile(hostFilename);
  string line;
  int tmpState, tmpVID, tmpInf, tmpHostK;
  vector<Host*> iniSusceptibles;
  vector<Host*> iniInfecteds;
  vector<Host*> iniRecovereds;
  map<int, Virus*> viruses = readViruses(virusFilename);
  Host* H;
  Virus* V;
  seedVirus = viruses[0];

  if(hostFile.is_open()){
    getline(hostFile,line);
    // Get each row from the host file
    while(getline(hostFile, line)){
      vector<int> tmpRow;
      stringstream ss(line);
      int i;
      while(ss >> i){
	tmpRow.push_back(i);
	if(ss.peek() == ',' || ss.peek() == '\n' || ss.peek() == ' ') ss.ignore();
      }
      tmpSize = tmpRow.size();
      if(tmpSize == 4){
	tmpState = tmpRow[0];
	tmpVID = tmpRow[1];
	tmpInf = tmpRow[2];
	tmpHostK = tmpRow[3];
	//	if(tmpState == 0) cout << tmpState << " " << tmpVID << " " << tmpInf << " " << tmpHostK << endl;
	H = new Host(static_cast<State>(tmpState), this, tmpHostK);
	// If this host was infected, get the infecting virus ID and get this virus from the viruses vector
	if(tmpState == Infected){
	  H->infect(viruses[tmpInf],viruses[tmpInf]->getBirth());
	  viruses[tmpInf]->updateHost(H);
	  iniInfecteds.push_back(H);
	}
	else if(tmpState == Susceptible){
	  if(tmpVID != 0){
	    H->addInfection(viruses[tmpVID]);
	    viruses[tmpVID]->updateHost(H);
	  }
	  iniSusceptibles.push_back(H);
	}
	else{
	  if(tmpVID != 0){
	    H->addInfection(viruses[tmpVID]);
	    viruses[tmpVID]->updateHost(H);
	  }
	  iniRecovereds.push_back(H);
	}
      }
    }
    hostFile.close();
  }
  susceptibles = iniSusceptibles;
  infecteds = iniInfecteds;
  recovereds = iniRecovereds;
  /*Rcpp::Rcout << "Virus ID: ";
  Virus::printIDgenerator();
  Rcpp::Rcout << endl;
  Rcpp::Rcout << "True Virus ID: " << viruses.rbegin()->first +1 << endl;*/
  Virus::set_generator(viruses.rbegin()->first+1);
  Rcpp::Rcout << "Read hosts complete" << endl;
  Rcpp::Rcout << "#########################" << endl;
}

std::map<int, Virus*> HostPopulation::readViruses(std::string virusFilename){
  Rcpp::Rcout << "-------------------------------------------------------" << endl;
  Rcpp::Rcout << "Reading in viruses... " << endl;
  std::map <int, Virus*> viruses;
  vector<int> parentIDs;
  ifstream virusFile(virusFilename);
  string line;
  Virus* V;
  int vid, birth, death, level, infectionNo;
  double bindingAvid, bindingAvidIni,distanceToParent, distanceToRoot, distHost, tmpK, changeV, changeR;
  if(virusFile.is_open()){
    getline(virusFile,line);
    while(getline(virusFile,line)){
      vector<double> tmpRow;
      stringstream ss(line);
      double i;
      while(ss >> i){
	tmpRow.push_back(i);
	if(ss.peek() == ',' || ss.peek() == '\n' || ss.peek() == ' ') ss.ignore();	
      }
      vid = (int)tmpRow[0];
      birth = (int)tmpRow[1];
      death = (int)tmpRow[2];
      parentIDs.push_back((int)tmpRow[3]);
      level = (int)tmpRow[4];
      bindingAvidIni = tmpRow[5];	
      bindingAvid = tmpRow[6];
      infectionNo = (int)tmpRow[7];
      distanceToParent = tmpRow[8];
      distanceToRoot = tmpRow[9];
      distHost = tmpRow[10];
      tmpK = tmpRow[11];      
      changeV = tmpRow[12];
      changeR = tmpRow[13];

      V = new Virus(vid, birth, death, NULL, level, bindingAvidIni, bindingAvid, infectionNo, distanceToParent, distanceToRoot, distHost, tmpK, NULL, changeV, changeR);
      viruses.insert(make_pair(vid,V));
    }
    typedef map<int, Virus*>::iterator virusIt;
    int index = 0;
    for(virusIt iterator = viruses.begin(); iterator != viruses.end(); iterator++){
      // Iterate through all viruses and update parents
      if(parentIDs[index] >= 0){
	iterator->second->updateParent(viruses[parentIDs[index]]);
      }
      index++;
    }
  }
  Rcpp::Rcout << "Read viruses complete" << endl;
  Rcpp::Rcout << "-------------------------------------------------------" << endl;
  return(viruses);
}

