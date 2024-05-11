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
    double c = 0.299792458;
    double tan_cerenkov = 0.900404; //42 degrees

    double t_c = par[0];
    double z_c = par[1];
    double coszenith = par[2];
    double d_c = par[3];
    
    return t_c + 1/c * (-(x[0]-z_c)*coszenith + tan_cerenkov*sqrt(d_c*d_c + (1-coszenith*coszenith)*(x[0]-z_c)*(x[0]-z_c)));
}

void fit::graph_fit(TTree* tree){
    
    std::vector<float> *t = new std::vector<float>();
    std::vector<float> *domid = new std::vector<float>();
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
        float min_t = *std::min_element(t->begin(), t->end());
        for (int j=0; j<n; j++){
            gr->SetPoint(j, 300-30*(domid->at(j)), t->at(j));
            // std::cout<<"z: "<<300-30*(domid->at(j))<<" t: "<<t->at(j)<<std::endl;
            zmin = std::min<double>(zmin, double(300-30*(domid->at(j))));
            zmax = std::max<double>(zmax, double(300-30*(domid->at(j))));
        }

        TF1 *f = new TF1("f", fit_function, zmin, zmax, 4);
        f->SetParameter(0, min_t);
        f->SetParameter(2, -0.5);
        gr->Fit(f, "");
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