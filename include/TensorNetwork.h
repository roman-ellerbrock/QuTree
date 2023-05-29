#pragma once
#include "Graph.h"
#include "backend/Tensor.h"

namespace qutree {

using NetworkShape = Graph<std::vector<index_t>>;

NetworkShape standardShape(const Graph<> &graph, index_t bondDimension,
                            index_t leafDimension);

/**
 * \brief Selects which part of a graph are active
 */
using GraphSelector = std::pair<bool, bool>;

/// @brief tensor network: has tensors at nodes but not at edges/leaves
extern GraphSelector TN;
/// @brief matrix network: has matrices at edges/leaves but not at nodes
extern GraphSelector MN;
/// @brief canonical tensor network: has tensors everywhere
extern GraphSelector CTN;

using TensorNetwork = Graph<Tensor>;

/// signature of functions that generate tensors like 
using Generator = tensorlib::Tensor (*)(tensorlib::IntArrayRef, tensorlib::TensorOptions);

template <class tn>
TensorNetwork createTN(const Graph<tn> &shape, GraphSelector s,
    Generator g = tensorlib::rand);

NetworkShape testTree();

} // namespace qutree


std::ostream& operator<<(std::ostream& os, const qutree::NetworkShape& shape);
std::ostream& operator<<(std::ostream& os, const qutree::TensorNetwork& tn);




