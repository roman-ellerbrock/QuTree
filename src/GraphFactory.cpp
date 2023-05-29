#include "GraphFactory.h"
#include "TensorNetwork.h"
#include <numeric>

namespace qutree {
bool contains(const std::vector<Edge> &edges, Edge edge) {
  for (Edge e : edges) {
    if (e == edge) {
      return true;
    }
  }
  return false;
}

template <class Attribute>
void addToGraph(Graph<Attribute> &subGraph, const Graph<Attribute> &graph,
                Edge edge) {
  std::vector<Edge> edges = graph.upEdges(edge);
  for (Edge e : edges) {
    subGraph.edges_[e] = graph.edges_.at(e);
    subGraph.edges_[flip(e)] = graph.edges_.at(flip(e));
    subGraph.nodes_[to(e)] = graph.edges_.at(e);
    addToGraph(subGraph, graph, e);
  }
}

/**
 * \brief Select smallest subgraph that connects all leaves
 **/
template <class Attribute>
Graph<Attribute> subgraph(const Graph<Attribute> &graph,
                          const std::vector<Node> &leaves) {
  Graph<Attribute> subgraph;
  for (Node leaf_id : leaves) {
    Edge leaf = graph.leafEdge(leaf_id);
    subgraph.edges_[leaf] = graph.edges_.at(leaf);
    subgraph.nodes_[to(leaf)] = graph.nodes_.at(to(leaf));
    addToGraph(subgraph, graph, leaf);
  }
  return subgraph;
}

template Graph<std::any> subgraph(const Graph<std::any> &graph,
                                  const std::vector<Node> &leaves);
template NetworkShape subgraph(const NetworkShape &graph,
                               const std::vector<Node> &leaves);

/**
 *        4
 *      /   \
 *    2      3
 *   / \
 *  0   1
 *
 */
/**
 *        6
 *      /   \
 *    4      5
 *   / \    / \
 *  0   1  2   3
 *
 */

void addLayer(Graph<> &graph, std::list<index_t> &buffer) {
  index_t idx = buffer.back() + 1;
  for (index_t i = 0; i < buffer.size(); i += 2) {
    // add a new node to graph and buffer
    graph.nodes_[idx] = "";
    buffer.push_back(idx);

    // Add first edge
    Edge edge = {buffer.front(),
                 idx}; // create edge and go to next list element
    graph.edges_[edge] = "";
    graph.edges_[flip(edge)] = ""; // inverse edge
    buffer.pop_front();

    // Add second edge
    edge = {buffer.front(), idx};
    graph.edges_[edge] = "";
    graph.edges_[flip(edge)] = ""; // inverse edge
    buffer.pop_front();
    idx++;
  }
}

template <class Attribute>
Graph<Attribute> balancedBinaryTree(index_t nLeaves) {
  Graph graph;

  std::list<index_t> buffer(nLeaves);
  std::iota(buffer.begin(), buffer.end(), 0);

  // add leaves
  for (auto i : buffer) {
    /// check!
    graph.edges_[Leaf({-(i + 1), i})] = "";
    graph.nodes_[Node(i)] = "";
  }

  while (buffer.size() > 1) {
    addLayer(graph, buffer);
  }
  return graph;
}
template Graph<std::any> balancedBinaryTree<std::any>(index_t nLeaves);

} // namespace qutree
