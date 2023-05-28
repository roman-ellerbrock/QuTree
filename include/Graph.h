#pragma once
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>

#include <iostream>
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
bool isInLeaf(Edge e);
bool isOutLeaf(Edge e);
bool isLeaf(Edge e);

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
 * Leaves are just edges that point from negative (logic) Nodes
 */
template <class Attribute = std::any> class Graph {
public:
  Graph() = default;
  ~Graph() = default;

  Graph(const std::map<Node, Attribute> &&nodes,
        const std::map<Edge, Attribute> &&edges)
      : nodes_(nodes), edges_(edges) {}

  Graph(const std::map<Node, Attribute> &nodes,
        const std::map<Edge, Attribute> &edges)
      : nodes_(nodes), edges_(edges) {}

  /// cross-initialization
  Graph(const Graph<Attribute> &graph)
      : nodes_(graph.nodes_), edges_(graph.edges_) {}

  template <class B> Graph(const Graph<B> &graph) {
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
  }

  bool empty() const {
    return (nodes_.empty() && edges_.empty());
  }

  std::vector<Edge> inEdges(Node node) const;

  std::vector<Edge> outEdges(Node node) const;

  std::vector<Edge> inEdges(Edge edge) const;

  std::vector<Edge> outEdges(Edge edge) const;

  std::vector<Edge> upEdges(Node node) const;

  std::vector<Edge> upEdges(Edge edge) const;

  std::vector<Edge> downEdges(Node node) const;

  std::vector<Edge> downEdges(Edge edge) const;

  /// @brief returns the index that an outgoing edge is attached to from(edge)
  /// @param node node that edge goes out from
  /// @param edge outgoing edge
  /// @return attachment index [0, nneighbors[
  index_t outIndex(Node node, Edge edge) const;

  /// @brief returns the index that an incoming edge is attachde to to(edge)
  /// @param node node that edge comes in to
  /// @param edge incoming edge
  /// @return attachment index [0, nneighbors[
  index_t inIndex(Edge edge, Node node) const;

  /**
   * \brief Generate vector of edges sorted so that it can be sweeped over in
   *TNS simulations
   **/
  std::vector<Edge> sortedEdges() const;

  bool containsNode(Node node) const {
    return (nodes_.find(node) != nodes_.end());
  }

  bool containsEdge(Edge edge) const {
    return (edges_.find(edge) != edges_.end());
  }

  Edge leafEdge(Node leaf) const {
    auto vec = outEdges(leaf);
    return vec.front();
  }

  Node root() const { return nodes_.rbegin()->first; }

  std::map<Node, Attribute> nodes_;
  std::map<Edge, Attribute> edges_;
};

index_t layer(Node node, const Graph<> &graph);

} // namespace qutree

std::ostream &operator<<(std::ostream &os, const qutree::Edge &edge);
std::ostream &operator<<(std::ostream &os, const qutree::Graph<> &graph);
