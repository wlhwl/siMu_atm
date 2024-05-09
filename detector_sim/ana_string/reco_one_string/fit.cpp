#include <math.h>
#include <algorithm>
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "fit.hpp"
#include "root_editor.hpp"


fit::fit(/* args */)
{
}

fit::~fit()
{
}

double fit::fit_function(double* x, double *par){
    double c = 299792458;
    double sin_cerenkov = 0.669; //42 degrees
    //par[0]=t0, par[1]=z_c, par[2]=sintheta = coszenith, par[3]=d_c, x=z
    return par[0] + ((x[0] - par[1]) * par[2] + sin_cerenkov * sqrt(pow(par[3],2) + (1 - pow(par[2],2)) * pow((x[0] - par[1]),2)))/c;
}

void fit::graph_fit(TTree* tree){
    
    std::vector<float> *t = new std::vector<float>();
    std::vector<int> *domid = new std::vector<int>();
    // std::vector<double> *zenith = new std::vector<double>();
    double coszenith;
    std::vector<float> *eventid = new std::vector<float>();

    TBranch* eventid_branch = root_editor::load_branch(tree, "EventId", &eventid);
    TBranch* t_branch = root_editor::load_branch(tree, "t", &t);
    TBranch* domid_branch = root_editor::load_branch(tree, "DomId", &domid);
    TBranch* zenith_branch = root_editor::load_branch(tree, "coszenith", &coszenith);


    int entries = tree->GetEntries();;
    for (int i=0; i<entries; i++){

        eventid_branch->GetEntry(i);
        std::cout<<"eventid: "<<eventid->at(0)<<std::endl;
        t_branch->GetEntry(i);
        domid_branch->GetEntry(i);
        zenith_branch->GetEntry(i);

        int n = domid->size();
        TGraph *gr = new TGraph(n);

        double zmin = -300;
        double zmax = 300;

        for (int j=0; j<n; j++){
            gr->SetPoint(j, 300-30*(domid->at(j)), t->at(j));
            zmin = std::min<double>(zmin, double(300-30*(domid->at(j))));
            zmax = std::max<double>(zmax, double(300-30*(domid->at(j))));
            
        }

        TF1 *f = new TF1("f", fit_function, zmin, zmax, 4);
        gr->Fit(f, "S");
        // zenith->push_back(f->GetParameter(2));
        coszenith = f->GetParameter(2);
        zenith_branch->Fill();
        delete gr;
        delete f;
    }
    tree->Write("", TObject::kOverwrite);
    delete t;
    delete domid;
    delete eventid;

}