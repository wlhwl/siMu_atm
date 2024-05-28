#include <math.h>
#include <algorithm>
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResult.h"
#include <map>
#include "../include/fit.hpp"
#include "../include/root_editor.hpp"


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
    double chi2;
    int ndf;
    double prob;
    int fitstatus;
    std::vector<float> *eventid = new std::vector<float>();

    TBranch* eventid_branch = root_editor::load_branch(tree, "EventId", &eventid);
    TBranch* t_branch = root_editor::load_branch(tree, "t", &t);
    TBranch* domid_branch = root_editor::load_branch(tree, "DomId", &domid);
    TBranch* zenith_branch = root_editor::load_branch(tree, "coszenith", &coszenith);
    TBranch* chi2_branch = root_editor::load_branch(tree, "chi2", &chi2);
    TBranch* ndf_branch = root_editor::load_branch(tree, "ndf", &ndf);
    TBranch* prob_branch = root_editor::load_branch(tree, "prob", &prob);
    TBranch* fitstatus_branch = root_editor::load_branch(tree, "fitstatus", &fitstatus);


    int entries = tree->GetEntries();;
    for (int i=0; i<entries; i++){

        eventid_branch->GetEntry(i);
        std::cout<<"eventid: "<<eventid->at(0)<<std::endl;
        t_branch->GetEntry(i);
        domid_branch->GetEntry(i);
        zenith_branch->GetEntry(i);

        std::map<float,std::vector<float>> domid_to_t;
        for (int j=0; j<domid->size(); j++){
            domid_to_t[domid->at(j)].push_back(t->at(j));
        }
        std::map<float,std::vector<float>>::iterator it = domid_to_t.begin();

        int n = domid_to_t.size();
        TGraphErrors* gr = new TGraphErrors(n);
        // TGraph *gr = new TGraph(n);

        double zmin = -300;
        double zmax = 300;
        double t_sigma = 1.3/(2*sqrt(2*log(2)));//ns  )
        float min_t = *std::min_element(t->begin(), t->end());
        
        for (int j=0; j<n; j++, it++){
            float DomId = it->first;
            float t_value = it->second[0];
            gr->SetPoint(j, 300-30*DomId, t_value);
            gr->SetPointError(j, 0, t_sigma);
            // std::cout<<"z: "<<300-30*(domid->at(j))<<" t: "<<t->at(j)<<std::endl;
            zmin = std::min<double>(zmin, double(300-30*DomId));
            zmax = std::max<double>(zmax, double(300-30*DomId));
        }

        TF1 *f = new TF1("f", fit_function, zmin, zmax, 4);
        f->SetParameter(0, min_t);
        f->SetParameter(2, 0.92);
        f->SetParLimits(2, 0, 1.0);
        //TODO check other parameters
        TFitResultPtr fitresult =  gr->Fit(f, "S");
        coszenith = fitresult->Parameter(2);
        chi2 = fitresult->Chi2();
        ndf = fitresult->Ndf();
        prob = fitresult->Prob();
        fitstatus = fitresult->Status();
        //TODO record uncertainty

        zenith_branch->Fill();
        chi2_branch->Fill();
        ndf_branch->Fill();
        prob_branch->Fill();
        fitstatus_branch->Fill();

        delete gr;
        delete f;
    }
    tree->Write("", TObject::kOverwrite);
    delete t;
    delete domid;
    delete eventid;
}