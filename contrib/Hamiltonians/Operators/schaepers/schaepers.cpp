//
// Created by hoppe on 04.03.21.
//

#include "schaepers.h"

namespace Operator{

    extern "C" {
        void hinitschaepers_(double* coeff, int* diag, int* modes, int* atomnum, double* masses, int* coup);
        void h_(int* mode, int* teil, double* hPsi, double* Psi,
                int* dim, double* matrix, double* trafo, double* ort);
    }

    schaepers::schaepers(const Tree& tree, vector<double> masses, vector<int> coupling) {

        if(masses.size() - 2 != coupling.size()){
            cerr << "ERROR: number of arguments supplied is not correct (schaepers KEO)" << endl;
            exit(1);
        }

        masses_ = masses;
        coupling_ = coupling;

        SpecialInitialize(tree);
    }

    void schaepers::callHinit(Vectorcd &coeffs, Matrix<int> &diag) {

        double* mass = new double[masses_.size()];
        for(int i = 0; i < masses_.size(); i++){
            mass[i] = masses_[i];
        }

        int* coup = new int[coupling_.size()];
        for(int i = 0; i < coupling_.size(); i++){
            coup[i] = coupling_[i];
        }

        int len = masses_.size() + 2;

        hinitschaepers_((double*)(&coeffs(0)),(int*)&diag(0,0),&nmodes,&len,mass,coup);

        delete[] mass;
        delete[] coup;
    }

    void schaepers::InitOperator() {
        nmodes = 12;
        nparts = 999999; // not needed really
        SystemH = h_;
    }
} // namespace Operators