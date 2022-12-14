//
// Created by hoppe on 3/18/21.
//

#include "liuch4cl.h"

extern "C" {
void ch4clpot_(double* v, double* q, double* m);
void potini_();
void pipnn_(double* v, double* q);
}


liuch4cl::liuch4cl(vector<double> massvec, vector<int> coupling) {
    potini_();
    massvec_ = massvec;
    coupling_ = coupling;
}

double liuch4cl::evaluate(const Vectord &Xv, size_t part) const {

    double v;
    std::array<double,12> q{};
    std::array<double,6> m{};

    m[0] = massvec_[1]; // C

    m[1] = massvec_[0]; // non-reactive H
    m[2] = massvec_[0];
    m[3] = massvec_[0];

    m[4] = massvec_[2]; // reactive H
    m[5] = massvec_[3]; // Cl

    for(int i = 0; i < 9; ++i){
        q[i] = Xv[i];
    }

    // move Cl
    for(int i = 9; i < 12; ++i){
        q[i] = 10000.;
    }


    ch4clpot_(&v,q.data(),m.data());
    const double vres = v;
//    std::cout << vres << std::endl;


    /*
    std::array<double,18> q;
    // currently: C,H,H,H,H_reac,Cl
    // target: H,H,H,H_reac,C,Cl

    // move CH3-Hs
    for(int i = 0; i < 9; ++i){
        q[i] = Xv[i + 3];
    }

    // place reactive H
    q[9] = 1000.;
    q[10] = 1000.;
    q[11] = 1000.;

    // move carbon
    q[12] = Xv[0];
    q[13] = Xv[1];
    q[14] = Xv[2];

    // place Cl
    q[15] = 1000. - 2.*1.27;
    q[16] = 1000.;
    q[17] = 1000.;

    double vres;
    for(auto& i : q){
        i *= 0.5291772;
    }
    pipnn_(q.data(),&vres);
    vres /= 27.2114;
    */


    return vres;
}
