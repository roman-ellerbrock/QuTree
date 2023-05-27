#pragma once
#include "Operator.h"

namespace qutree {

namespace pauli {
Tensor X();
Tensor Y();
Tensor Z();
} // namespace pauli

SumOfProducts transversalFieldIsing(double J, double g, index_t N);

} // namespace qutree