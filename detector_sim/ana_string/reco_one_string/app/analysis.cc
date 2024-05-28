#include <fstream>
#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "../include/fit.hpp"
#include "../include/filter.hpp"
#include "../include/process_g4_npe.hpp"
#include "../include/root_editor.hpp"

int main(int argc, char** argv){
    
    std::shared_ptr<TFile> root_file = std::make_shared<TFile>("copy.root", "UPDATE");
    if(root_file->IsZombie()){
        std::cerr << "Error: cannot open file " << "copy.root" << std::endl;
        return -1;
    }

    TTree* tree = root_editor::load_tree(root_file.get(),"string_tree");
 
    // process_g4_npe npe_process(tree);
    // npe_process.set_qe_file("../qe.csv");
    // npe_process.set_tts(1.3);
    // npe_process.process();
    
    filter filter_obj(tree);
    filter_obj.apply_L2(5);
    TTree* filtered_tree = filter_obj.get_filtered_tree();

    fit::graph_fit(filtered_tree);
   
    root_file->Close();
    return 0;
}                                                                           