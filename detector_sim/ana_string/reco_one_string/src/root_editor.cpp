#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "trident/root_editor.hpp"

root_editor::root_editor() {
    // Constructor
}

root_editor::~root_editor() {
    // Destructor
}

TTree* root_editor::load_tree(TFile* file, std::string treename){
    TTree* tree = dynamic_cast<TTree *>(file->Get(treename.c_str()));
    if(tree == nullptr){
        std::cerr << "Tree " << treename << " not found in file" << std::endl;
    }
    return tree;
}
