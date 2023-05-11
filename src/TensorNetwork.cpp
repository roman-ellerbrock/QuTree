#include "TensorNetwork.h"

namespace qutree {

GraphSelector TN = {true, false, false};
GraphSelector MT = {false, true, true};

NetworkShape standardShape(const Graph<> &graph, index_t bondDimension,
                           index_t leafDimension) {
  NetworkShape shape(graph);

  for (auto &nodepair : shape.nodes_) {
    auto &ref = nodepair.second;
    Node node = nodepair.first;
    auto nNeighbors = graph.outEdges(node).size();
    ref = std::vector<index_t>(nNeighbors, bondDimension);
  }

  for (auto &edgepair : shape.edges_) {
    auto &ref = edgepair.second;
    ref = {bondDimension, bondDimension};
  }

  for (auto &leafpair : shape.leaves_) {
    auto &ref = leafpair.second;
    ref = {leafDimension, leafDimension};
  }

  return shape;
}

TensorNetwork randomTTNS(const NetworkShape &shape, GraphSelector s) {
  TensorNetwork psi(shape);
  using namespace std;
  bool nodes = get<0>(s);
  bool edges = get<1>(s);
  bool leaves = get<2>(s);

  if (nodes) {
    for (auto &nodepair : psi.nodes_) {
      Node node = nodepair.first;
      auto &ref = nodepair.second;
      ref = tensorlib::rand(shape.nodes_.at(node), tensorlib::options());
    }
  }

  if (edges) {
    for (auto &edgepair : psi.edges_) {
      Edge edge = edgepair.first;
      auto &ref = edgepair.second;
      ref = tensorlib::rand(shape.edges_.at(edge), tensorlib::options());
    }
  }

  if (leaves) {
    for (auto &leafpair : psi.leaves_) {
      Leaf leaf = leafpair.first;
      auto &ref = leafpair.second;
      ref = tensorlib::rand(shape.leaves_.at(leaf), tensorlib::options());
    }
  }

  return psi;
}

TensorNetwork createTensorNetwork(const Graph<> &graph) {

  //  TensorNetwork psi(graph);
  TensorNetwork psi;

  return psi;
}

} // namespace qutree