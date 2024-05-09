#include <fstream>
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include <random>
#include "../filter.hpp"
#include "../root_editor.hpp"

int main(int argc, char** argv){
    
    std::shared_ptr<TFile> root_file = std::make_shared<TFile>("copy.root", "UPDATE");
    if(root_file->IsZombie()){
        std::cerr << "Error: cannot open file " << "copy.root" << std::endl;
        return -1;
    }

    TTree* tree = root_editor::load_tree(root_file.get(),"string_tree");

    filter filter_obj(tree);
    filter_obj.apply_L2(4);
    
    root_file->Close();
    
    return 0;
}                             