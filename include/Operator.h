#pragma once
#include "Graph.h"
#include "TensorNetwork.h"
#include "backend/Tensor.h"

namespace qutree {

using ProductOperator = std::map<index_t, Tensor>;
using SumOfProducts = std::vector<std::pair<scalar_t, ProductOperator>>;

std::vector<Node> leaves(const ProductOperator& P);

TensorNetwork matrixNetwork(const NetworkShape &shape, const ProductOperator &P);

} // namespace qutree