#pragma once
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>

#include "yaml-cpp/yaml.h"
#include <iostream>
#include <nlohmann/json.hpp>

#include "backend/Tensor.h"
#include <any>

namespace qutree {

using Node = index_t;               // uid
using Leaf = std::pair<Node, Node>; // uid
using Edge = std::pair<Node, Node>; // {from, to}

Node from(Edge e);
Node to(Edge e);
Edge flip(Edge e);
bool isUpEdge(Edge e);
bool isDownEdge(Edge e);

template <class Attribute>
std::vector<Edge> mapToVector(const std::map<Edge, Attribute> &m);
void eraseValue(std::vector<Edge> &edges, Edge e);

/**
 * \class Graph
 * \brief This class represents a graph consisting of nodes, edges, and leaves.
 *
 * @param Attribute Defines what is associates with nodes, edges, and leaves.
 *
 * The graph class encodes the topology of a tensor network. Nodes, edges and
 * leaves are stored in maps. They can be given names which are stored in the
 * maps as well.
 */
template <class Attribute = std::any> class Graph {
public:
  Graph() = default;
  ~Graph() = default;

  Graph(const std::map<Node, Attribute> &&nodes,
        const std::map<Edge, Attribute> &&edges,
        const std::map<Leaf, Attribute> &&leaves)
      : nodes_(nodes), edges_(edges), leaves_(leaves) {}

  Graph(const std::map<Node, Attribute> &nodes,
        const std::map<Edge, Attribute> &edges,
        const std::map<Leaf, Attribute> &leaves)
      : nodes_(nodes), edges_(edges), leaves_(leaves) {}

  /// cross-initialization
  Graph(const Graph<Attribute> &graph)
      : nodes_(graph.nodes_), edges_(graph.edges_), leaves_(graph.leaves_) {}

  template <class B> Graph(const Graph<B> &graph) {
    /// Leaves
    for (const auto &p : graph.leaves_) {
      Leaf l = p.first;
      leaves_[l] = Attribute();
    }

    /// Edges
    for (const auto &p : graph.edges_) {
      Edge e = p.first;
      edges_[e] = Attribute();
    }

    /// Nodes
    for (const auto &p : graph.nodes_) {
      Node n = p.first;
      nodes_[n] = Attribute();
    }
  }

  void clear() {
    nodes_.clear();
    edges_.clear();
    leaves_.clear();
  }

  bool empty() const {
    return (nodes_.empty() && edges_.empty() && leaves_.empty());
  }

  std::vector<Edge> inEdges(Node node) const {
    std::map<Edge, Attribute> filtered_edges;
    std::copy_if(edges_.begin(), edges_.end(),
                 std::inserter(filtered_edges, filtered_edges.end()),
                 [node](const std::pair<Edge, Attribute> &p) {
                   return to(p.first) == node;
                 });

    return mapToVector(filtered_edges);
  }

  std::vector<Edge> outEdges(Node node) const {
    std::map<Edge, Attribute> filtered_edges;
    std::copy_if(edges_.begin(), edges_.end(),
                 std::inserter(filtered_edges, filtered_edges.end()),
                 [node](const std::pair<Edge, Attribute> &p) {
                   return (from(p.first) == node);
                 });

    return mapToVector(filtered_edges);
  }

  std::vector<Edge> inEdges(Edge edge) const {
    std::vector<Edge> prev = inEdges(from(edge));
    eraseValue(prev, flip(edge));
    return prev;
  }

  std::vector<Edge> outEdges(Edge edge) const {
    std::vector<Edge> post = outEdges(to(edge));
    eraseValue(post, flip(edge));
    return post;
  }

  std::vector<Edge> upEdges(Node node) const {
    std::map<Edge, Attribute> filtered_edges;
    std::copy_if(edges_.begin(), edges_.end(),
                 std::inserter(filtered_edges, filtered_edges.end()),
                 [node](const std::pair<Edge, Attribute> &p) {
                   return (from(p.first) == node) && isUpEdge(p.first);
                 });

    return mapToVector(filtered_edges);
  }

  std::vector<Edge> upEdges(Edge edge) const { return upEdges(to(edge)); }

  std::vector<Edge> downEdges(Node node) const {
    std::map<Edge, Attribute> filtered_edges;
    std::copy_if(edges_.begin(), edges_.end(),
                 std::inserter(filtered_edges, filtered_edges.end()),
                 [node](const std::pair<Edge, Attribute> &p) {
                   return (from(p.first) == node) && isDownEdge(p.first);
                 });

    return mapToVector(filtered_edges);
  }

  std::vector<Edge> downEdges(Edge edge) const { return downEdges(to(edge)); }

  /**
   * \brief Generate vector of edges sorted so that it can be sweeped over in
   *TNS simulations
   **/
  std::vector<Edge> sortedEdges() const {
    std::map<Edge, Attribute> filtered_edges;
    /// copy up-edges, i.e. where from(edge) < to(edge)
    std::copy_if(edges_.begin(), edges_.end(),
                 std::inserter(filtered_edges, filtered_edges.end()),
                 [](const std::pair<Edge, Attribute> &p) {
                   return (isUpEdge(p.first));
                 });
    std::vector<Edge> sorted_edges = mapToVector(filtered_edges);

    /// copy down-edges, i.e. where from(edge) > to(edge)
    filtered_edges.clear();
    std::copy_if(edges_.begin(), edges_.end(),
                 std::inserter(filtered_edges, filtered_edges.end()),
                 [](const std::pair<Edge, Attribute> &p) {
                   return (isDownEdge(p.first));
                 });
    std::vector<Edge> down = mapToVector(filtered_edges);
    std::reverse(down.begin(), down.end());
    sorted_edges.insert(sorted_edges.end(), down.begin(), down.end());
    return sorted_edges;
  }

  bool containsNode(Node node) const {
    return (nodes_.find(node) != nodes_.end());
  }

  bool containsEdge(Edge edge) const {
    return (edges_.find(edge) != edges_.end());
  }

  bool containsLeaf(Leaf leaf) const {
    return (leaves_.find(leaf) != leaves_.end());
  }

  Node root() const { return nodes_.rbegin()->first; }

  std::map<Node, Attribute> nodes_;
  std::map<Edge, Attribute> edges_;
  std::map<Leaf, Attribute> leaves_;
};

index_t layer(Node node, const Graph<> &graph);

} // namespace qutree

std::ostream &operator<<(std::ostream &os, const qutree::Edge &edge);
std::ostream &operator<<(std::ostream &os, const qutree::Graph<> &graph);