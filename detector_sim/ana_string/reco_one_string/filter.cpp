#include "TTree.h"
#include "filter.hpp"
#include "root_editor.hpp"
#include <iostream>
#include <set>

filter::filter(TTree* t){
    tree = t;
    tree_br.branches[0] = root_editor::load_branch(tree,"EventId",&tree_vec.EventId);
    tree_br.branches[1] = root_editor::load_branch(tree,"PmtId",&tree_vec.PmtId);
    tree_br.branches[2] = root_editor::load_branch(tree,"DomId",&tree_vec.DomId);
    tree_br.branches[3] = root_editor::load_branch(tree,"pe",&tree_vec.pe);
    tree_br.branches[4] = root_editor::load_branch(tree,"t0",&tree_vec.t0);
    tree_br.branches[5] = root_editor::load_branch(tree,"e0",&tree_vec.e0);
    tree_br.branches[6] = root_editor::load_branch(tree,"wavelength",&tree_vec.wavelength);
    tree_br.branches[7] = root_editor::load_branch(tree,"t",&tree_vec.t);

}

filter::~filter(){
}

void filter::create_filtered_tree(std::string& tree_name){
    //create a filtered tree
    filtered_tree = new TTree(tree_name.c_str(),tree_name.c_str());
    filtered_tree->SetDirectory(tree->GetDirectory());

    fil_br.branches[0] = root_editor::load_branch(filtered_tree,"EventId",&fil_vec.EventId);
    fil_br.branches[1] = root_editor::load_branch(filtered_tree,"PmtId",&fil_vec.PmtId);
    fil_br.branches[2] = root_editor::load_branch(filtered_tree,"DomId",&fil_vec.DomId);
    fil_br.branches[3] = root_editor::load_branch(filtered_tree,"wavelength",&fil_vec.wavelength);
    fil_br.branches[4] = root_editor::load_branch(filtered_tree,"t",&fil_vec.t);

}

void filter::apply_L2(int dom_num_thres){
    std::string filtered_tree_name = "L2_" + std::to_string(dom_num_thres) + "doms";
    create_filtered_tree(filtered_tree_name);
    
    for(int i = 0; i < tree->GetEntries(); i++){
        root_editor::branches_getentry(tree_br,i);
        root_editor::branches_getentry(fil_br,i);
        std::set<float> doms;
        for(int j = 0; j < tree_vec.DomId->size(); j++){
            if(tree_vec.pe->at(j) > 0){
                doms.insert(tree_vec.DomId->at(j));
                fil_vec.EventId->push_back(tree_vec.EventId->at(j));
                fil_vec.PmtId->push_back(tree_vec.PmtId->at(j));
                fil_vec.DomId->push_back(tree_vec.DomId->at(j));
                fil_vec.wavelength->push_back(tree_vec.wavelength->at(j));
                fil_vec.t->push_back(tree_vec.t->at(j));
            }
        }
        if(doms.size() >= dom_num_thres){
            filtered_tree->Fill();
        }
        fil_vec.EventId->clear();
        fil_vec.PmtId->clear();
        fil_vec.DomId->clear();
        fil_vec.wavelength->clear();
        fil_vec.t->clear();
    }
    filtered_tree->Write("",TObject::kOverwrite);
}

void filter::first_hit(){
    
}