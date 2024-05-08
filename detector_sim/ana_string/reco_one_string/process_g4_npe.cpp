#include <fstream>
#include <iostream>
#include <sstream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include <random>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include "process_g4_npe.hpp"
#include "root_editor.hpp"

process_g4_npe::process_g4_npe(TTree* t){
    tree = t;
}

process_g4_npe::~process_g4_npe(){

}

void process_g4_npe::process(){

    if(qe_file != "")
        apply_qe();
    
    if(tts != 0)
        apply_tts(); 
}

void process_g4_npe::apply_qe() {
    //load branches
    std::vector<float> *wavelength = new std::vector<float>;
    TBranch* wavelength_branch = root_editor::load_branch(tree,"wavelength",&wavelength);
    std::vector<float> *e0 = new std::vector<float>;
    TBranch* e0_branch = root_editor::load_branch(tree,"e0",&e0);
    std::vector<int> *pe = new std::vector<int>;
    TBranch* pe_branch = root_editor::load_branch(tree,"pe",&pe);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 1);

    QEvalue qe = get_qe(qe_file);

    int entries = tree->GetEntries();

    tree->SetAutoSave(100);

    for (int i=0; i<entries; i++){
        pe->clear();
        wavelength->clear();

        wavelength_branch->GetEntry(i);
        e0_branch->GetEntry(i);
        pe_branch->GetEntry(i);
        
        for (int j=0;j<e0->size();j++){
            wavelength->push_back(1240/(e0->at(j)));
            if (dis(gen) < (interpolate(qe,wavelength->at(j))/100.0)){
                pe->push_back(1);
            }else{
                pe->push_back(0);
            }
        }
        wavelength_branch->Fill();
        pe_branch->Fill();
    
    }
    tree->AutoSave();
    tree->Write("",TObject::kOverwrite);
    
    delete e0;
    delete wavelength;
    delete pe;
}

void process_g4_npe::apply_tts() {
    //load branches
    std::vector<float> *t = new std::vector<float>;
    TBranch* t_branch = root_editor::load_branch(tree, "t", &t);
    std::vector<float> *t0 = new std::vector<float>;
    TBranch* t0_branch = root_editor::load_branch(tree, "t0", &t0);
    std::vector<int> *coolpe = new std::vector<int>;
    TBranch* pe_branch = root_editor::load_branch(tree, "pe", &coolpe);
    
    int entries = tree->GetEntries();
    
    tree->SetAutoSave(100);

    float sigma = tts / (2*sqrt(2*log(2)));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<float> dis(0, sigma);
    
    for (int i=0; i<entries; i++){
        t->clear();

        t_branch->GetEntry(i);
        t0_branch->GetEntry(i);
        pe_branch->GetEntry(i);

        for (int j=0; j <t0->size(); j++){
            if (coolpe->at(j) == 1){
                t->push_back(t0->at(j) + dis(gen));
            }else{
                t->push_back(0);
            }
        }
        t_branch->Fill();
    }
    tree->AutoSave();
    tree->Write("",TObject::kOverwrite);
    
    delete t;
    delete t0;
    delete coolpe;
}

float process_g4_npe::interpolate(QEvalue QE, float wl){
    // TODO: check if the interpolator is right
    boost::math::interpolators::cardinal_cubic_b_spline<double> spline(QE.qe.begin(), QE.qe.end(), QE.wl.front(),QE.wl[1]-QE.wl[0]);
    return spline(wl);
}

QEvalue process_g4_npe::get_qe(std::string& qe_file){
    
    QEvalue qe;
    std::ifstream file(qe_file);
    if (!file.is_open()){
        std::cerr << "Could not open file " << qe_file << std::endl;
        return qe;
    }
    
    std::string line;
    //negelect the first row and check it btw
    if(!std::getline(file,line)){
        std::cerr << "Could not read file " << qe_file << std::endl;
        return qe;
    }

    while (std::getline(file,line)){
        std::istringstream ss(line);
        std::string field;
        //wl
        if(!std::getline(ss, field, ',')){
            std::cerr << "Could not read wavelength from file " << qe_file << std::endl;
            return qe;
        }
        qe.wl.push_back(std::stof(field));
        //qe
        if(!std::getline(ss, field, ',')){
            std::cerr << "Could not read qe from file " << qe_file << std::endl;
            return qe;
        }
        qe.qe.push_back(std::stof(field));
    }
    file.close();
    return qe;
}
