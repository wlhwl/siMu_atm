#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include <iostream>
#include <fstream>
#include <vector>


template <typename T,typename U>
void load_branch(TTree *oldtree, std::string oldbranchname, T oldaddress, TTree *newtree, std::string newbranchname, U newaddress){
    TBranch *branch = oldtree->GetBranch(oldbranchname.c_str());
    if(branch==nullptr){
        std::cerr << "Error: cannot find branch " << oldbranchname << std::endl;
        return;
    }
    branch->SetAddress(oldaddress);
    // std::cout<<typeid(address).name()<<std::endl;
    if(!newtree->FindBranch(newbranchname.c_str())){
        newtree->Branch(newbranchname.c_str(), newaddress);
    }
}

int main(int argc, char **argv){
    if (argc !=2){
        std::cerr << "Usage: " << argv[0] << " <String id>" << std::endl;
        return 1;
    }
    //create the string root file
    std::string string_file_path = "string_" + std::string(argv[1]) + ".root";
    std::shared_ptr<TFile> string_file = std::make_shared<TFile>(string_file_path.c_str(), "RECREATE");
    if (string_file->IsZombie()){
        std::cerr << "Error: cannot open file " << string_file_path << std::endl;
        return -1;
    }
    //create the string tree
    std::shared_ptr<TTree> string_tree = std::make_shared<TTree>("string_tree", ("string id " + std::string(argv[1])).c_str());
    string_tree->SetDirectory(string_file.get());
    //list branches name
    std::string branch_name[] = {"t0","e0","PmtId","DomId"};
    std::vector<float> *data_t0, *data_e0;
    std::vector<int> *data_pmtid, *data_domid;
    std::vector<float> *string_t0, *string_e0, *string_pmtid, *string_domid;
    data_t0 = new std::vector<float>;
    data_e0 = new std::vector<float>;
    data_pmtid = new std::vector<int>;
    data_domid = new std::vector<int>;
    string_t0 = new std::vector<float>;
    string_e0 = new std::vector<float>;
    string_pmtid = new std::vector<float>;
    string_domid = new std::vector<float>;

    //record event id
    std::vector<float> *string_eventid = new std::vector<float>;
    string_tree->Branch("EventId", &string_eventid);

    for(int jobid=0;jobid<24;jobid++){

        std::string job_file_path = "../../..//single_muon_jobs/job_" + std::to_string(jobid) + "/copy.root";
        std::shared_ptr<TFile> job_file = std::make_shared<TFile>(job_file_path.c_str(), "READ");
        if (job_file->IsZombie()){
            std::cerr << "Error: cannot open file " << job_file_path << std::endl;
            return -1;
        }
        TTree *job_tree = dynamic_cast<TTree *>(job_file->Get("PmtHit"));
        if (job_tree == nullptr){
            std::cerr << "Error: cannot find tree PmtHit in file " << job_file_path << std::endl;
            return -1;
        }
        //load branches
        load_branch(job_tree, branch_name[0], &data_t0, string_tree.get(), branch_name[0], &string_t0);
        load_branch(job_tree, branch_name[1], &data_e0, string_tree.get(), branch_name[1], &string_e0);
        load_branch(job_tree, branch_name[2], &data_pmtid, string_tree.get(), branch_name[2], &string_pmtid);
        load_branch(job_tree, branch_name[3], &data_domid, string_tree.get(), branch_name[3], &string_domid);

        int entries = job_tree->GetEntries();

        for(int ientry=0;ientry<entries;ientry++){
            data_t0->clear();
            data_e0->clear();
            data_pmtid->clear();
            data_domid->clear();
            string_t0->clear();
            string_e0->clear();
            string_pmtid->clear();
            string_domid->clear();
            string_eventid->clear();

            job_tree->GetEntry(ientry);
            int vector_size = data_domid->size();
            // std::cout<<vector_size<<std::endl;
            if(vector_size==0){
                continue;
            }
            for(int j=0;j<vector_size;j++){
                if(int((data_domid->at(j))/21)==std::stoi(argv[1])){
                    string_t0->push_back(data_t0->at(j));
                    string_e0->push_back(data_e0->at(j));
                    string_pmtid->push_back(static_cast<float>(data_pmtid->at(j)));
                    string_domid->push_back(static_cast<float>(data_domid->at(j)));
                    string_eventid->push_back(static_cast<float>(ientry+jobid*10260));
                }
            }
            if(string_domid->size()!=0){
                string_tree->Fill();
            }else{
                continue;
            }
        }
        string_file.get()->cd();
        string_tree->Write("",TObject::kOverwrite);
        job_file->Close();
    }
    string_file->Close();
    return 0;
}
