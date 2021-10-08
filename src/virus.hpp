#ifndef VIRUS_HPP
#define VIRUS_HPP

#include "hostpopulation.hpp"

class Host;

class Virus {
private:
    static int v_IDgenerator; // For indexing new viruses
    static int n_variants; // How many variants are there?
    static Rcpp::IntegerVector variantGenerators; // Track how many of each variant we've made
    
    static Rcpp::NumericMatrix crossImmunity; // Matrix of cross-immunity between viruses
    static Rcpp::NumericMatrix vlPars; // Store a pre-computed matrix of randomly generate viral kinetics parameters. Pre-computed for speed. 
    static Rcpp::NumericMatrix infectiousnessPars; // Matrix with a row for each variant, giving parameters of the VL/infectiousness relationship
    // Need to be careful with the size, as we don't want to create more viruses then there are entries. Otherwise loop round to the start.
    
    static double true_0; // True 0 viral load value
    
    int id; // Unique ID of this virus
    int variant_id; // Unique ID for this variant
    int variant; // Variant type, used for working out cross-reactivity etc
    
    int birth; // Time of creation
    int death; // Time of recovery
    int infectionNo; // From perspective of the host
    
    Virus* parent; // Pointer to parent virus
    Host* host; // Pointer to host
    
    // Viral kinetics parameters of this infection
    double tg;
    double tp;
    double to;
    double tw;
    double alpha;
    double symptomatic;
    
    // Infectiousness parameters
    double infectiousnessMax;
    double infectiousnessGradient;
    double infectiousnessInflection;
    
public:
    // Constructors
    Virus();
    
    // Constructor for simulation
    Virus(int _t, int _variant, Virus* _parent, Host* _host, int _infectionNo);
    
    // Destructor
    ~Virus(){}; 

    // ID generator functions
    static void printIDgenerator();
    static int getIDgenerator();
    static void set_generator(int _start);
    static void initiateVariantGenerator();
    static int getVariantGenerator(int variant);
    static void setVariantGenerator(int variant, int _start);
    
    // Class-wide parameters
    static void set_vlPars(Rcpp::NumericMatrix new_vlPars);
    static void set_infectiousnessPars(Rcpp::NumericMatrix new_infectiousnessPars);
    static void set_crossImmunity(Rcpp::NumericMatrix new_crossImmunity);
    static void set_variants(int _n_variants);
    static void set_true_0(double _true_0);
    
    // Class-wide functions
    static double getCrossImmunity(Virus* infectingVirus, Virus* immuneVirus);
    
    // Infection state
    void updateParent(Virus* newParent);
    void updateHost(Host* newHost);
    void kill(int cur_t);
    
    // Set viral kinetics parameters
    void set_default();
    void set_tg(double new_tg);
    void set_tp(double new_tp);
    void set_to(double new_to);
    void set_tw(double new_tw);
    void set_alpha(double new_alpha);
    void set_variant(int new_variant);
    void set_symptomatic(bool new_symptomatic);
    
    // Set infectiousness parameters
    void set_infectiousnessMax(double new_set_infectiousnessMax);
    void set_infectiousnessGradient(double new_set_infectiousnessGradient);
    void set_infectiousnessInflection(double new_infectiousnessInflection);
 
    // Access states and parameters   
    int getId();
    int getVariantId();
    int getBirth();
    int getDeath();
    int getInfectionNo();
    int getVariant();
    int getAgeOfParentAtBirth();
    
    bool hasRecovered(int cur_t);
    bool hasOnset(int cur_t);
    
    Host* getHost();
    Virus* getParent(); 
    
    // Access kinetics parameters
    double get_tg();
    double get_tp();
    double get_to();
    double get_tw();
    double get_alpha();
    bool get_symptomatic();
    
    // Access infectiousness parameters
    double get_infectiousnessMax();
    double get_infectiousnessGradient();
    double get_infectiousnessInflection();
    
    // Calculations/events
    double calculateViralLoad(int cur_t);
    double getInfectiousness(int cur_t);
};

#endif
