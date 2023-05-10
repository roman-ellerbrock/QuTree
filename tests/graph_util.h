#pragma once
#include "Graph.h"

namespace qutree {

std::map<Node, std::any> create_nodes();

std::map<Edge, std::any> create_edges();

std::map<Leaf, std::any> create_leaves();

Graph<> binary_4_graph();

} // namespace qutree