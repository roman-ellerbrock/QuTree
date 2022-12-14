//
// Created by Roman Ellerbrock on 4/28/22.
//

#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H
#include "TreeOperators/Potential.h"
#include "TreeOperators/CoordinateTransformation.h"
#include "TreeOperators/Hamiltonian.h"

Vectord getx0(const Tree& tree, const CoordinateTransformation& U);

void find_minimum(const Hamiltonian& H, const Tree& tree);

#endif //NORMAL_MODES_H
