#pragma once
#include "Graph.h"
#include "backend/Tensor.h"

namespace qutree {

using ProductOperator = std::map<index_t, Tensor>;
using SumOfProducts = std::vector<std::pair<scalar_t, ProductOperator>>;

} // namespace qutree