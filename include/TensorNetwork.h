#pragma once
#include "Graph.h"
#include "backend/Tensor.h"

namespace qutree {

/**
 *
 */

using NetworkShape = Graph<std::vector<index_t>>;

NetworkShape standardShape(const Graph<> &graph, index_t bondDimension,
                            index_t leafDimension);

/**
 * \brief Selects which part of a graph are active
 */
using GraphSelector = std::tuple<bool, bool, bool>;

/// @brief tensor network: has tensors at nodes but not at edges/leaves
extern GraphSelector TN;
/// @brief matrix network: has matrices at edges/leaves but not at nodes
extern GraphSelector MT;

using TensorNetwork = Graph<Tensor>;

TensorNetwork randomTTNS(const NetworkShape &shape, GraphSelector s);

TensorNetwork createTensorNetwork(const Graph<> &graph);

} // namespace qutree
