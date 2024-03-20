#ifndef INITIAL_ENERGY_GENERATOR
#define INITIAL_ENERGY_GENERATOR

#include "TRandom3.h"
#include "TF1.h"


class initial_energy_generator{
public:
    initial_energy_generator(double e_min, double e_max, double power_index, unsigned int r_seed)
        : energy_min(e_min), energy_max(e_max), seed(r_seed), power_index(power_index){
        rnd = new TRandom3(seed);
        gRandom->SetSeed(seed);
        f_classical_logE = new TF1("classical logE",classical_logE_distribution, log(energy_min), log(energy_max), 1);
        f_classical_logE->SetParameter(0, power_index);
        f_logE_minus1 = new TF1("(logE)^-1","1/x",log(energy_min), log(energy_max));
    }

    ~initial_energy_generator(){
        delete rnd;
        delete f_classical_logE;
    }

    double get_E_classical_distribution(){
        double logE = f_classical_logE->GetRandom(log(energy_min), log(energy_max));
        return exp(logE);
    }

    double get_E_expUniform(){
        double power = rnd->Uniform(log10(energy_min), log10(energy_max));
        return pow(10,power);
    }

    double get_E_minus1(){
        double logE = f_logE_minus1->GetRandom(log(energy_min), log(energy_max));
        return exp(logE);
    }

private:
    double energy_min;
    double energy_max;
    unsigned int seed;
    double power_index;
    TRandom3* rnd;
    TF1* f_classical_logE;
    TF1* f_logE_minus1;

    static Double_t classical_logE_distribution(Double_t *x, Double_t *par){
        // // energy distribution: pow(E,-2.7) for E<3e15eV; pow(E,-3.1) otherwise
        // // f(E)=pow(E,a) ~ f(lnE)=exp( lnE * (a+1) )
        // // log(3e15) = 35.637389; exp(-1.7*35.637389)/exp(-2.1*35.637389) = 1551845.6
        // if(*x<35.637389) return exp((*x)*(-1.7));
        // else return 1551845.6*exp((*x)*(-2.1));
        return exp((*x)*(par[0]+1.0));
    }
};


#endif
