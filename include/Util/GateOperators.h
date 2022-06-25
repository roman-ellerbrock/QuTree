//
// Created by Roman Ellerbrock on 8/27/21.
//

#ifndef GATEOPERATORS_H
#define GATEOPERATORS_H
#include "Core/Matrix.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeOperators/SOPVector.h"

/**
 * Rationale:
 * - Small library for Quantum Circuits
 * - Only contains functions that are needed for unit testing
 */

Matrixcd sigma_x();
Matrixcd sigma_y();
Matrixcd sigma_z();

/// Turn an operation into a controlled operation
SOPcd controlled(const SOPcd& in, size_t control);
SOPcd controlled(const MLOcd& M, size_t c, size_t t);
SOPcd controlled(const Matrixcd& L, size_t c, size_t t);

SOPcd CNot(size_t c, size_t t);
SOPcd CZ(size_t c, size_t t);

SOPVectorcd QFT(size_t mode_start, size_t n_bit, bool adjungate = false, size_t approx = 0);

#endif //GATEOPERATORS_H
