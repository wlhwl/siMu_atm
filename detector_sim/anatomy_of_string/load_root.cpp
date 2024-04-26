#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <fstream>

int main(int argc, char **argv){

    if(argc!=2){
        std::cerr << "Usage: " << argv[0] << " <ROOT_file>" << std::endl;
        return -1;
    }

    std::string ROOTPath = "../single_muon_jobs/job_" + std::string(argv[1]) + "/copy.root";
    TFile *file = new TFile(ROOTPath.c_str());
    if(file->IsZombie()){
        std::cerr << "Error: cannot open file " << ROOTPath << std::endl;
        return -1;
    }
    TTree *tree = (TTree*)file->Get("DomHit");
    if(tree==nullptr){
        std::cerr << "Error: cannot find tree DomHit in file " << ROOTPath << std::endl;
        return -1;
    }
    std::ofstream csv_survived,csv_dead;
    csv_survived.open("survive_event.csv", std::ios_base::app);
    csv_dead.open("dead_event.csv", std::ios_base::app);
        std::vector<float> *id=0;
        tree->SetBranchAddress("DomId", &id);
        for(int entry=0; entry<tree->GetEntries(); entry++){
            tree->GetEntry(entry);
            if(id->size()==0){
                csv_dead << 10260*std::stoi(argv[1])+entry << std::endl;
            }else{
                for(int i : *id){
                    csv_survived << 10260*std::stoi(argv[1])+entry << "," << i << std::endl;
                    }
            }
        }
        csv_dead.close();
        csv_survived.close();
        file->Close();
        return 0;
    }
