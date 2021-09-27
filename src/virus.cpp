#include "virus.hpp"
#include "host.hpp"
#include<iostream>
#include<math.h>
#include<vector>
#include <Rcpp.h>

using namespace std;

int Virus::v_IDgenerator = 1;
double Virus::true_0 = 0;

Virus::Virus(int _id, int _birth, int _death, int _variant, Virus* _parent, Host* _host, int _infectionNo,
             double _tg, double _tp, double _to, double _tw, double _alpha
             ){
    id = _id;
    birth = _birth;
    death = _death;
    variant = _variant;
    infectionNo = _infectionNo;
    
    parent = _parent;
    host = _host;
    
    tg = _tg;
    tp = _tp;
    to = _to;
    tw = _tw;
    alpha = _alpha;
}

Virus::Virus(int _t, int _variant, Virus* _parent, Host* _host, int _infectionNo,
             double _tg, double _tp, double _to, double _tw, double _alpha
){
    id = v_IDgenerator++;
    birth = _t;
    death = -1;
    variant = _variant;
    infectionNo = _infectionNo;
    
    parent = _parent;
    host = _host;
    
    tg = _tg;
    tp = _tp;
    to = _to;
    tw = _tw;
    alpha = _alpha;
}


Virus::Virus(int _id, int _birth, int _death, int _variant, Virus* _parent, Host* _host){
    id = _id;
    birth = _birth;
    death = _death;
    variant = _variant;
    
    parent = _parent;
    host = _host;
    
    tg = 2;
    tp = 3;
    to = 5;
    tw = 10;
    alpha = 10;
}

void Virus::set_default(){
    tg = 2;
    tp = 3;
    to = 5;
    tw = 10;
    alpha = 10;
}

// ID generator functions
void Virus::printIDgenerator(){
    Rcpp::Rcout << v_IDgenerator;
}

int Virus::getIDgenerator(){
    return v_IDgenerator;
}

void Virus::set_generator(int _start){
    v_IDgenerator = _start;
}

// Status and property access functions
int Virus::getId(){
    return id;
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

int Virus::getInfectionNo(){
    return infectionNo;
}

// Set states and parameters
void Virus::kill(int cur_t){
    death = cur_t;
}

void Virus::updateParent(Virus* newParent){
    parent = newParent;
}

void Virus::updateHost(Host* newHost){
    host = newHost;
}

void Virus::set_tg(double new_tg){tg=new_tg;}
void Virus::set_tw(double new_tw){tw=new_tw;}
void Virus::set_tp(double new_tp){tp=new_tp;}
void Virus::set_to(double new_to){to=new_to;}
void Virus::set_alpha(double new_alpha){alpha=new_alpha;}
void Virus::set_variant(int new_variant){variant=new_variant;}


// Viral load calculation
double calculateViralLoad(int cur_t){
    int t = cur_t - birth;
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
        y = t*alpha/tp;
    } else if(t > tp){
        y = alpha - (t-tp)*(alpha/(to - tp + tw));
    } else {
        y = true_0;
    }
    return y;
    
}


