//
// Created by hoppe on 3/18/21.
//

#ifndef MCTDH_CPP_LIUCH4CL_H
#define MCTDH_CPP_LIUCH4CL_H

#pragma once
#include "TreeOperators/Potential.h"

class liuch4cl : public Potential {
public:
    liuch4cl() = default;
    liuch4cl(vector<double> massvec, vector<int> coupling);

    ~liuch4cl() = default;

    double evaluate(const Vectord& Xv, size_t part) const override;

private:
    vector<double> massvec_;
    vector<int> coupling_;

};


#endif //MCTDH_CPP_LIUCH4CL_H
