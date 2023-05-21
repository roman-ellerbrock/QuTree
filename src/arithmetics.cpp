#include "arithmetics.h"
#include <exception>

namespace qutree {

TensorNetwork dotProduct(const TensorNetwork &Bra, const TensorNetwork &Ket) {
  //    TensorNetwork mt = createTN(Bra, MN, tensorlib::eyeWrapper);
  TensorNetwork mt = createTN(Bra, MN, tensorlib::eyeWrapper);

  for (Edge edge : mt.sortedEdges()) {
    if (isLeaf(edge)) {
      continue;
    }
    Node node = from(edge);
    Tensor &mat = mt.edges_[edge];
    Tensor bra = Bra.nodes_.at(node);
    Tensor ket = Ket.nodes_.at(node);

    // multiply in every previous edge's matrix
    std::vector<Edge> inEdges = mt.inEdges(edge);
    for (auto e : inEdges) {
      if (isLeaf(e)) {
        continue;
      }
      index_t idx = mt.inIndex(e, node);
      ket = tensorlib::contractMatrixTensor(mt.edges_[e], ket, idx);
    }
    // contract the tensors into the new edge-matrix
    index_t idx = mt.outIndex(node, edge);
    mat = tensorlib::contractTensorTensor(bra, ket, idx);
  }

  return mt;
}
} // namespace qutree
