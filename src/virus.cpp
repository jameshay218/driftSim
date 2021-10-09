#include "virus.hpp"
#include "host.hpp"
#include<iostream>
#include<math.h>
#include<vector>
#include <Rcpp.h>

using namespace std;

/* ========================================== */
// Class-specific variables
/* ========================================== */
int Virus::v_IDgenerator = 1; // For indexing new viruses
int Virus::n_variants = 1; // How many variants exist?
double Virus::true_0 = 0; // True 0 viral load value

Rcpp::IntegerVector Virus::variantGenerators; // Track how many of each variant we've made

Rcpp::NumericMatrix Virus::crossImmunity ; // Matrix of cross-immunity between viruses
Rcpp::NumericMatrix Virus::vlPars; // Matrix of all vl pars to draw from
Rcpp::NumericMatrix Virus::infectiousnessPars; // Matrix of infectiousness parameters, row for each variant

/* ========================================== */
// Constructors
/* ========================================== */
Virus::Virus(int _t, int _variant, Virus* _parent, Host* _host, int _infectionNo){
    id = v_IDgenerator++; // Increment overall virus generator
    variant_id = variantGenerators(_variant)++; // Increment generator for this variant
    
    birth = _t;
    death = -1;
    variant = _variant;
    infectionNo = _infectionNo;
    
    parent = _parent;
    host = _host;
    
    int index;
    
    if(vlPars.nrow() > 0){
        // A little funky but feels efficient to code. We have a pre-computed matrix of values to pull from.
        // We have a matrix for each variant. There are the same number of rows for each variant. These
        // matrices get concatenated and we pull from the correct "chunk" of the combined matrix.
        // A vector of indexers keeps track of how many of each variant we've made.
        index = variant*(vlPars.nrow()/n_variants) + variant_id % (vlPars.nrow()/n_variants);
        tg = vlPars(index, 0);
        tp = vlPars(index, 1);
        to = vlPars(index, 2);
        tw = vlPars(index, 3);
        alpha = vlPars(index, 4);
        symptomatic = vlPars(index, 5);
        
    } else {
        tg = 2;
        tp = 3;
        to = 5;
        tw = 10;
        alpha = 10;
        symptomatic=false;
    }
    
    if(infectiousnessPars.nrow() > 0){
        infectiousnessMax = infectiousnessPars(variant, 0);
        infectiousnessGradient = infectiousnessPars(variant, 1);
        infectiousnessInflection = infectiousnessPars(variant, 2);
    } else {
        infectiousnessMax = 0.5;
        infectiousnessGradient = 0.75;
        infectiousnessInflection = 5;
    }
}

Virus::Virus(double _t, double _tg, double _tp, double _to, double _tw, double _alpha, 
             double _infectiousnessMax, double _infectiousnessGradient, double _infectiousnessInflection){
    birth = _t;
    death = -1;
    variant = 0;
    infectionNo = 0;
    
    parent = NULL;
    host = NULL;
    
    tg = _tg;
    tp = _tp;
    to = _to;
    tw = _tw;
    alpha = _alpha;
    symptomatic=false;
    
    infectiousnessMax = _infectiousnessMax;
    infectiousnessGradient = _infectiousnessGradient;
    infectiousnessInflection = _infectiousnessInflection;
}

// Default constructor
Virus::Virus(){
    id = v_IDgenerator++;
    variant_id = variantGenerators(variant)++; // Increment generator for this variant
    
    birth = 0;
    death = -1;
    variant = 0;
    infectionNo = 0;
    
    parent = NULL;
    host = NULL;
    
    tg = 2;
    tp = 3;
    to = 5;
    tw = 10;
    alpha = 10;
    symptomatic=false;
    
    infectiousnessMax = 0.5;
    infectiousnessGradient = 0.75;
    infectiousnessInflection = 5;
}


/* ========================================== */
// ID generator functions
/* ========================================== */
// Overall ID
void Virus::printIDgenerator(){
    Rcpp::Rcout << v_IDgenerator;
}

int Virus::getIDgenerator(){
    return v_IDgenerator;
}

void Virus::set_generator(int _start){
    v_IDgenerator = _start;
}

// Variant-specific IDs
void Virus::initiateVariantGenerator(){
    variantGenerators = Rcpp::NumericVector(n_variants);
}

int Virus::getVariantGenerator(int variant){
    return(variantGenerators(variant));
}

void Virus::setVariantGenerator(int variant, int _start){
    variantGenerators(variant) = _start;
}

/* ========================================== */
// Manage class-wide parameters
/* ========================================== */
void Virus::set_vlPars(Rcpp::NumericMatrix new_vlPars){
    vlPars = new_vlPars;
}

void Virus::set_infectiousnessPars(Rcpp::NumericMatrix new_infectiousnessPars){
    infectiousnessPars = new_infectiousnessPars;
}

void Virus::set_crossImmunity(Rcpp::NumericMatrix new_crossImmunity){
    crossImmunity = new_crossImmunity;
}

void Virus::set_variants(int _n_variants){
    n_variants = _n_variants;
}
void Virus::set_true_0(double _true_0){
    true_0 = _true_0;
}

/* ========================================== */
// Class-wide functions
/* ========================================== */
// Find cross-immunity between two viruses based on the Virus::crossImmunity matrix
double Virus::getCrossImmunity(Virus* infectingVirus, Virus* immuneVirus){
    return crossImmunity(infectingVirus->getVariant(), immuneVirus->getVariant());
}

/* ========================================== */
// Set individual parameters
/* ========================================== */
void Virus::updateParent(Virus* newParent){
    parent = newParent;
}

void Virus::updateHost(Host* newHost){
    host = newHost;
}

void Virus::kill(double cur_t){
    death = cur_t;
}

// Set viral kinetics parameters
void Virus::set_default(){
    tg = 2.0;
    tp = 3.0;
    to = 5.0;
    tw = 10.0;
    alpha = 10.0;
    symptomatic=false;
    
    infectiousnessMax = 0.5;
    infectiousnessGradient = 0.75;
    infectiousnessInflection = 5;
}

void Virus::set_tg(double new_tg){tg=new_tg;}
void Virus::set_tw(double new_tw){tw=new_tw;}
void Virus::set_tp(double new_tp){tp=new_tp;}
void Virus::set_to(double new_to){to=new_to;}
void Virus::set_alpha(double new_alpha){alpha=new_alpha;}
void Virus::set_symptomatic(bool new_symptomatic){symptomatic=new_symptomatic;}
void Virus::set_variant(int new_variant){variant=new_variant;}

// Set infectiousness parameters
void Virus::set_infectiousnessMax(double new_set_infectiousnessMax){infectiousnessMax=new_set_infectiousnessMax;}
void Virus::set_infectiousnessGradient(double new_set_infectiousnessGradient){infectiousnessGradient=new_set_infectiousnessGradient;}
void Virus::set_infectiousnessInflection(double new_infectiousnessInflection){infectiousnessInflection=new_infectiousnessInflection;}

/* ========================================== */
// Access states and parameters
/* ========================================== */
int Virus::getId(){
    return id;
}
int Virus::getVariantId(){
    return variant_id;
}
int Virus::getBirth(){
    return birth;
}

int Virus::getDeath(){
    return death;
}

Host* Virus::getHost(){
    return host;
}
Virus* Virus::getParent(){
    return parent;
}

int Virus::getInfectionNo(){
    return infectionNo;
}

int Virus::getVariant(){
    return variant;
}

int Virus::getAgeOfParentAtBirth(){
    if(parent->getBirth() >= 0){
        return birth - parent->getBirth();
    }
    return(-1);
}

double Virus::get_tg(){return tg;}
double Virus::get_tw(){return tw;}
double Virus::get_tp(){return tp;}
double Virus::get_to(){return to;}
double Virus::get_alpha(){return alpha;}
bool Virus::get_symptomatic(){return symptomatic;}

// Set infectiousness parameters
double Virus::get_infectiousnessMax(){return infectiousnessMax;}
double Virus::get_infectiousnessGradient(){return infectiousnessGradient;}
double Virus::get_infectiousnessInflection(){return infectiousnessInflection;}


// If current time is long enough after creation to have 0 viral load, then should be recovered
bool Virus::hasRecovered(double cur_t){
    if(cur_t - birth >= tp + tw){
        return true;
    } else {
        return false;
    }
}

// Check if onset is today
bool Virus::hasOnset(double cur_t){
    if(symptomatic &&  cur_t - birth == to){
        return true;
    }
    return false;
}

/* ========================================== */
// Model functions
/* ========================================== */
// Viral load calculation
double Virus::calculateViralLoad(double cur_t){
    double t = cur_t - birth;
    double y;
    
    // If virus has been killed
    if(cur_t > death && death > -1){
        y = true_0;
        return y;
    }
    
    // Pre-infection and pre-latent
    if(t <= tg){
        y = true_0;
    } else if(t > tg & t <= tp) {
        //y = (t-tg)*alpha/(tp-tg);
        y = (t-tg)*alpha/(tp-tg);
    } else if(t > tp){
        y = alpha - (t-tp)*(alpha/(to - tp + tw));
    } else {
        y = true_0;
    }
    
    // Can't have negative viral load
    y = std::max(y, true_0);
    
    return y;
}

// Get infectiousness right now based on current viral load and some assumed relationship between viral load
// and probability of transmitting
double Virus::getInfectiousness(double cur_t){
    double vl = calculateViralLoad(cur_t);
    
    // Just one assumption. Assume logistic curve for relationship between viral load and probability of transmission
    double infectiousness = infectiousnessMax*(1.0 - 1.0/(1.0 + exp(infectiousnessGradient*(vl - infectiousnessInflection))));
    return(infectiousness);
}

