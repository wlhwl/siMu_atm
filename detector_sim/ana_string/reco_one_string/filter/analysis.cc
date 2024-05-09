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

    std::cout<<"Filtering tree..."<<std::endl;
    filter filter_obj(tree);
    std::cout<<"class created"<<std::endl;
    filter_obj.apply_L2(4);
    std::cout<<"L2 applied"<<std::endl;
    root_file->Close();
    
    return 0;
}                             