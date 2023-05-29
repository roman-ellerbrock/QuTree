#include "Operator.h"
#include "GraphFactory.h"

namespace qutree {

TensorNetwork matrixNetwork(const NetworkShape &graph,
                            const ProductOperator &P) {

  std::vector<Node> ls = leaves(P);
  auto sgraph = subgraph(graph, ls);

  TensorNetwork mn = createTN(sgraph, MN);
  for (const auto &pairs : P) {
    index_t m = pairs.first;
    Tensor X = pairs.second;
    Edge leaf = graph.leafEdge(m);
    mn.edges_[leaf] = X;
  }
  return mn;
}

std::vector<Node> leaves(const ProductOperator &P) {
  std::vector<Node> ls;
  for (auto x : P) {
    ls.push_back(x.first);
  }
  return ls;
}

} // namespace qutree