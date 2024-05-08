#include "TTree.h"

class fit{
    public:
        fit(/* args */);
        ~fit();
        static double fit_function(double* x, double *par);
        static void graph_fit(TTree* );
};
