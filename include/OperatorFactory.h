#pragma once
#include "Operator.h"

namespace qutree {

namespace pauli {
Tensor X();
Tensor Y();
Tensor Z();
} // namespace pauli

Tensor identity(index_t n);
ProductOperator identity(const NetworkShape& shape);

SumOfProducts transversalFieldIsing(double J, double g, index_t N);

} // namespace qutree