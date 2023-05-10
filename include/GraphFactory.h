#include "Graph.h"
#include "backend/Tensor.h"

namespace qutree {

template <class Attribute>
Graph<Attribute> subgraph(const Graph<Attribute> &graph,
                          const std::vector<Leaf> &leaves);

template <class Attribute = std::any>
Graph<Attribute> balancedBinaryTree(index_t nLeaves);

} // namespace qutree
