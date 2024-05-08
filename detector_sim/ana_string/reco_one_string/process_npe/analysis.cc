#include <fstream>
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include <random>
#include "../process_g4_npe.hpp"
#include "../root_editor.hpp"

int main(int argc, char** argv){
    
    std::shared_ptr<TFile> root_file = std::make_shared<TFile>("copy.root", "UPDATE");
    if(root_file->IsZombie()){
        std::cerr << "Error: cannot open file " << "copy.root" << std::endl;
        return -1;
    }

    TTree* tree = root_editor::load_tree(root_file.get(),"string_tree");

    process_g4_npe npe_process(tree);
    
    npe_process.set_qe_file("../qe.csv");
    npe_process.set_tts(1.3);
    npe_process.process();

    root_file->Close();
    return 0;
}                                                                           