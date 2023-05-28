#include "Operator.h"
#include "GraphFactory.h"

namespace qutree {

template <class attr>
TensorNetwork tensorNetwork(const Graph<attr> &graph, const ProductOperator &P,
                            GraphSelector s) {
  /**
   * todo(Roman) subgraph uses leaves (edges) which doesn't make sense in the
   * current context todo(Roman) make subgraph use vector<index_t>
   * {-1, -3, -5} -> subgraph -> mats -> represent
   */
//  auto sgraph = subgraph(graph, leaves(P));
  auto mats = createTN(graph, MN);
  for (const auto &pairs : P) {
    index_t m = pairs.first;
    Tensor X = pairs.second;
    Edge leaf = graph.leafEdge(m);
    mats.edges_[leaf] = X;
  }
  return mats;
}

template TensorNetwork tensorNetwork(const Graph<Tensor> &shape,
                                     const ProductOperator &P, GraphSelector s);
template TensorNetwork tensorNetwork(const Graph<std::vector<index_t>> &shape,
                                     const ProductOperator &P, GraphSelector s);

std::vector<index_t> leaves(const ProductOperator &P) {
  std::vector<index_t> ls;
  for (auto x : P) {
    ls.push_back(x.first);
  }
  return ls;
}

} // namespace qutree