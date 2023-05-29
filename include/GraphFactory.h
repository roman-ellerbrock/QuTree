#pragma once
#include "Graph.h"
#include "backend/Tensor.h"

namespace qutree {

// new definitions

// create a subgraph connecting different leaves
template <class Attribute>
Graph<Attribute> subgraph(const Graph<Attribute> &graph,
                          const std::vector<Node> &leaves);

// create a (close-to) balanced binary tree
template <class Attribute = std::any>
Graph<Attribute> balancedBinaryTree(index_t nLeaves);

} // namespace qutree
