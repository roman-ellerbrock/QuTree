//
// Created by Roman Ellerbrock on 2/2/20.
//

#ifndef RANDOMMATRICES_H
#define RANDOMMATRICES_H
#include "Core/Matrix.h"
#include <random>

namespace RandomMatrices {
	Matrixcd GUE(size_t dim, mt19937& gen);
	Matrixd GOE(size_t dim, mt19937& gen);
}

#endif //RANDOMMATRICES_H
