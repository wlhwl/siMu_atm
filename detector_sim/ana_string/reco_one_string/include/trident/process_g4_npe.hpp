#include "TTree.h"

struct QEvalue{
    std::vector<float> wl;
    std::vector<float> qe;
};

class process_g4_npe{
    public:
        process_g4_npe(TTree* );
        ~process_g4_npe();

        void set_qe_file(std::string f){qe_file=f;};
        void set_tts(double t){tts=t;};

        void process();

        void apply_qe();
        void apply_tts();
        static float interpolate(QEvalue , float );
        static QEvalue get_qe(std::string& );

    private:
        TTree* tree;
        std::string qe_file = "";
        double tts=0;
};
