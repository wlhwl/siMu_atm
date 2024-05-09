#include "TTree.h"
#include <iostream>

class filter {
public:
    filter(TTree* );
    ~filter();
    void create_filtered_tree(std::string& );
    void apply_L2(int );

    TTree* get_filtered_tree(){
        return filtered_tree;
    };
    
private:
    TTree* tree;
    TTree* filtered_tree;
    
    struct filtered_values{
        std::vector<float>* EventId = nullptr;
        std::vector<float>* PmtId = nullptr;
        std::vector<float>* DomId = nullptr;
        std::vector<float>* wavelength = nullptr;
        std::vector<float>* t = nullptr;
    };
    struct filtered_branches{
        std::vector<TBranch*> branches;
        filtered_branches() : branches(5, nullptr) {}
    };
    struct tree_values{
        std::vector<float>* EventId = nullptr;
        std::vector<float>* PmtId = nullptr;
        std::vector<float>* DomId = nullptr;
        std::vector<float>* pe = nullptr;
        std::vector<float>* t0 = nullptr;
        std::vector<float>* e0 = nullptr;
        std::vector<float>* wavelength = nullptr;
        std::vector<float>* t = nullptr;
    };
    struct tree_branches{
        std::vector<TBranch*> branches;
        tree_branches() : branches(8, nullptr) {}
    };
    
    filtered_values fil_vec;
    filtered_branches fil_br;
    tree_values tree_vec;
    tree_branches tree_br;
};