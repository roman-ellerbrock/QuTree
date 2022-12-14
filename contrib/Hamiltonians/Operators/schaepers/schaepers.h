//
// Created by hoppe on 04.03.21.
//
// This is an implementation of the kinetic energy operator derived by Schaepers, Zhao and Manthe (https://doi.org/10.1016/j.chemphys.2018.02.025)
// to study reactions involving methane molecules
//
// to build this operator, some input is needed:
// vector<int> coupling:    determines the type of coupling for additional atoms used (stereographic = 0, spherical = 1)
// vector<double> mass:     is the vector of masses for the atoms in the order: 3 non-reactive H, reactive H, methyl-C, coupled atoms
// bool methylfixed:
// bool frozen:

// the order of coordinates is given as:
// rho, phi_rho, theta_rho, theta, phi, chi, r, s, t, R, S, T, ....
// 0    1        2          3      4    5    6  7  8  9  10 11 ...

#pragma once
#include "TreeOperators/SumOfProductsOperator.h"
#include <vector>
#include "TreeOperators/FortranSOP.h"


namespace Operator {

    //SOPcd schaepers(const vector<int>& coupling, const vector<double>& mass, bool methylfixed, bool frozen);

    class schaepers : public FortranSOP {
    public:
        schaepers(const Tree& tree, vector<double> masses, vector<int> coupling);
        ~schaepers() = default;

    private:
        void callHinit(Vectorcd& coeffs, Matrix<int>& diag) override;

        void InitOperator() override;

    vector<double> masses_;
    vector<int> coupling_;



    };


} // namespace operator

