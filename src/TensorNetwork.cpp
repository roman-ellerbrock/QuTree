#include "TensorNetwork.h"

namespace qutree {

GraphSelector TN = {true, false};
GraphSelector MN = {false, true};
GraphSelector CTN = {true, true};

NetworkShape standardShape(const Graph<> &graph, index_t bondDimension,
                           index_t leafDimension) {
  NetworkShape shape(graph);

  for (auto &nodepair : shape.nodes_) {
    auto &ref = nodepair.second;
    Node node = nodepair.first;
    auto edges = graph.inEdges(node);
    ref.resize(edges.size());
    for (auto e : edges) {
      index_t idx = shape.inIndex(e, node);
      if (isInLeaf(e)) {
        ref[idx] = leafDimension;
      } else if (isOutLeaf(e)) {
        ref[idx] = 1;
      } else {
        ref[idx] = bondDimension;
      }
    }
  }

  for (auto &edgepair : shape.edges_) {
    Edge e = edgepair.first;
    auto &ref = edgepair.second;
    if (isInLeaf(e)) {
      ref = {leafDimension, leafDimension};
    } else if (isOutLeaf(e)) {
      ref = {1, 1};
    } else {
      ref = {bondDimension, bondDimension};
    }
  }

  return shape;
}

const tensorlib::IntArrayRef sizes(const tensorlib::IntArrayRef &A) {
  return A;
}

const tensorlib::IntArrayRef sizes(const tensorlib::Tensor &A) {
  return A.sizes();
}

template <class tn>
TensorNetwork
createTN(const Graph<tn> &shape, GraphSelector s,
         tensorlib::Tensor (*function)(tensorlib::IntArrayRef,
                                       tensorlib::TensorOptions)) {
  TensorNetwork psi(shape);
  using namespace std;
  bool nodes = get<0>(s);
  bool edges = get<1>(s);

  if (nodes) {
    for (auto &nodepair : psi.nodes_) {
      Node node = nodepair.first;
      auto &ref = nodepair.second;
      ref = function(sizes(shape.nodes_.at(node)), tensorlib::options());
    }
  }

  if (edges) {
    for (auto &edgepair : psi.edges_) {
      Edge edge = edgepair.first;
      auto &ref = edgepair.second;
      ref = function(sizes(shape.edges_.at(edge)), tensorlib::options());
    }
  }

  return psi;
}

template TensorNetwork
createTN(const NetworkShape &shape, GraphSelector s,
         tensorlib::Tensor (*function)(tensorlib::IntArrayRef,
                                       tensorlib::TensorOptions));

template TensorNetwork
createTN(const TensorNetwork &A, GraphSelector s,
         tensorlib::Tensor (*function)(tensorlib::IntArrayRef,
                                       tensorlib::TensorOptions));

} // namespace qutree

std::ostream& operator<<(std::ostream& os, const qutree::NetworkShape& graph) {
  using namespace qutree;
  os << "Nodes:\n";
  for (auto p : graph.nodes_) {
    Node node = p.first;
    auto A = p.second;
    os << node << "\n" << A << std::endl;
  }

  os << "Edges:\n";
  for (Edge e : graph.sortedEdges()) {
    auto A = graph.edges_.at(e);
    os << e << "\n" << A << std::endl;
  }

  return os;
}

std::ostream& operator<<(std::ostream& os, const qutree::TensorNetwork& tn) {
  using namespace qutree;
  os << "Nodes:\n";
  for (auto p : tn.nodes_) {
    Node node = p.first;
    Tensor A = p.second;
    os << "address: " << node << " | " << A.sizes() << "\n" << A << std::endl;
  }

  os << "Edges:\n";
  for (Edge e : tn.sortedEdges()) {
    Tensor A = tn.edges_.at(e);
    os << "address: " << e << " | " << A.sizes() << "\n" << A << std::endl;
  }

  return os;
}
