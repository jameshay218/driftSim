#ifndef VIRUS_HPP
#define VIRUS_HPP

#include "hostpopulation.hpp"

class Host;

class Virus {
private:
    static int v_IDgenerator;
    static double true_0;
    
    int id; // Unique ID of this virus
    
    int variant; // Variant ID, used for working out cross-reactivity etc
    
    int birth; // Time of creation
    int death; // Time of recovery
    int infectionNo; // From perspective of the host
    
    Virus* parent;
    Host* host;
    
    // Viral kinetics parameters of this infection
    double tg;
    double tp;
    double to;
    double tw;
    double alpha;
    
    
public:
    static double getCrossImmunity(Virus* A, Virus* B);
    
    // Constructors
    Virus();
    
    // Constructor for simulation
    Virus(int _t, int _variant, Virus* _parent, Host* _host, int _infectionNo,
          double _tg, double _tp, double _to, double _tw, double _alpha);
    
    // Full constructor
    Virus(int _id, int _birth, int _death, int _variant, Virus* _parent, Host* _host,
          double _tg, double _tp, double _to, double _tw, double _alpha);
    
    // Constructor determinstic kinetics
    Virus(int _id, int _birth, int _death, int _variant, Virus* _parent, Host* _host);
    
    ~Virus(){}; 
    
    // Accessing attributes
    int getId();
    int getBirth();
    int getDeath();
    int getInfectionNo();
    int getVariant();
    
    Host* getHost();
    
    Virus* getParent(); 
    
    void updateParent(Virus* newParent);
    void updateHost(Host* newHost);
    
    // Calculations/events
    double calculateViralLoad(int cur_t);
    void kill(int cur_t);
    
    // Change static member variables
    static void printIDgenerator();
    static int getIDgenerator();
    static void set_generator(int _start);
    static void set_default();
    static void set_tg(double new_tg);
    static void set_tp(double new_tp);
    static void set_to(double new_to);
    static void set_tw(double new_tw);
    static void set_alpha(double new_alpha);
    static void set_variant(int new_variant);
};

#endif
