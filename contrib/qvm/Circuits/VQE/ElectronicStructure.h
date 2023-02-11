//
// Created by Roman Ellerbrock on 7/27/21.
//

#ifndef ELECTRONICSTRUCTURE_H
#define ELECTRONICSTRUCTURE_H
#include "TreeOperators/SumOfProductsOperator.h"
#include "JordanWigner.h"

JordanWigner::TwoIndex readTwoIndexIntegral(ifstream& file, size_t nOrbital);
JordanWigner::FourIndex readFourIndexIntegral(ifstream& file, size_t nOrbital);

//JordanWigner::TwoIndex readTwoIndexIntegral(const string& filename);
//JordanWigner::FourIndex readFourIndexIntegral(const string& filename);

Matrixd convertTwoIndex(const JordanWigner::TwoIndex& h);
Tensord convertFourIndex(const JordanWigner::FourIndex& h);

SOPcd electronicStructure(const string& filename);

#endif //ELECTRONICHAMILTONIAN_H
