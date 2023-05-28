#include "arithmetics.h"
#include <exception>

using namespace std;

namespace qutree {

void dotProduct(TensorNetwork &mt, const TensorNetwork &Bra,
                const TensorNetwork &Ket) {

  for (Edge edge : mt.sortedEdges()) {
    if (isLeaf(edge)) {
      continue;
    }
    /**
     * ! The problem is that outIdx gives 0 although inIdx is 0 for a different
     * ! edge. We have to only use in- or out-index or fix the issue that
     * in-/out- ! index are compatible. ? For now, I will use in-index on
     * flipped edge
     *
     * ? I think it makes sense to change how outIdx works: it should actually
     * use ? the in-edges and use their index. todo(Roman): change that.
     */
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
    index_t idx = mt.inIndex(flip(edge), node); // dadum
    mat = tensorlib::contractTensorTensor(bra, ket, idx);
  }
}
} // namespace qutree
