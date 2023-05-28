#pragma once
#include "Graph.h"
#include "TensorNetwork.h"
#include "backend/Tensor.h"

namespace qutree {

using ProductOperator = std::map<index_t, Tensor>;
using SumOfProducts = std::vector<std::pair<scalar_t, ProductOperator>>;

std::vector<index_t> leaves(const ProductOperator& P);

template <class attr>
TensorNetwork tensorNetwork(const Graph<attr> &shape, const ProductOperator &P,
                            GraphSelector s);

} // namespace qutree