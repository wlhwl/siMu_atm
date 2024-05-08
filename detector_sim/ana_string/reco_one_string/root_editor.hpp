#ifndef ROOT_EDITOR_HPP
#define ROOT_EDITOR_HPP 1

#include "TBranch.h"
#include "TFile.h"

class root_editor {
public:
    root_editor();
    ~root_editor();
    static TTree* load_tree(TFile* , std::string );
    
    template <typename T>
    static TBranch* load_branch(TTree* tree, std::string branch_name, T branch_data){
        TBranch* branch = nullptr;
        if(tree == nullptr){
            std::cerr << "Tree is null" << std::endl;
            return nullptr;
        }
        if(tree->FindBranch(branch_name.c_str())){
            branch = tree->GetBranch(branch_name.c_str());
            if (branch == nullptr) {
                std::cerr << "Branch " << branch_name << " not found in tree" << std::endl;
                return nullptr;
            }
            branch->SetAddress(branch_data);
        }else{
            branch = tree->Branch(branch_name.c_str(), branch_data);
        }
        return branch;
    }


};

#endif
