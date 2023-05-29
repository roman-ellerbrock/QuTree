#include "Graph.h"
#include "backend/Tensor.h"

namespace qutree {

Node from(Edge e) { return e.first; }
Node to(Edge e) { return e.second; }
Edge flip(Edge e) { return {to(e), from(e)}; }
bool isInLeaf(Edge e) { return (from(e) < 0); }
bool isOutLeaf(Edge e) { return (to(e) < 0); }
bool isLeaf(Edge e) { return isInLeaf(e) || isOutLeaf(e); };

std::ostream &operator<<(std::ostream &os, const Edge &e) {
  os << "{" << e.first << ", " << e.second << "}";
  return os;
}

template <class Attribute>
std::vector<Edge> mapToVector(const std::map<Edge, Attribute> &m) {
  std::vector<Edge> v;
  std::transform(m.begin(), m.end(), std::back_inserter(v),
                 [](const std::pair<Edge, Attribute> &p) { return p.first; });
  return v;
}

template std::vector<Edge>
mapToVector(const std::map<Edge, std::vector<index_t>> &m);
template std::vector<Edge>
mapToVector(const std::map<Edge, tensorlib::Tensor> &m);

template <class Attribute>
std::vector<Edge> Graph<Attribute>::inEdges(Node node) const {
  std::map<Edge, Attribute> filtered_edges;
  // copy every edge for that to(edge) == node
  std::copy_if(edges_.begin(), edges_.end(),
               std::inserter(filtered_edges, filtered_edges.end()),
               [node](const std::pair<Edge, Attribute> &p) {
                 return to(p.first) == node;
               });

  return mapToVector(filtered_edges);
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::outEdges(Node node) const {
  std::map<Edge, Attribute> filtered_edges;
  // copy every edge for that from(edge) == node
  std::copy_if(edges_.begin(), edges_.end(),
               std::inserter(filtered_edges, filtered_edges.end()),
               [node](const std::pair<Edge, Attribute> &p) {
                 return (from(p.first) == node);
               });

  return mapToVector(filtered_edges);
}

void eraseValue(std::vector<Edge> &edges, Edge e) {
  while (true) {
    auto it = std::find(edges.begin(), edges.end(), e);
    if (it == edges.end()) {
      return;
    }
    edges.erase(it);
  }
}

bool isUpEdge(Edge e) { return (to(e) > from(e)); }

bool isDownEdge(Edge e) { return !isUpEdge(e); }

index_t layer(Node node, const Graph<> &graph) {
  index_t l = 0;
  auto ps = graph.upEdges(node);
  while (!ps.empty()) {
    if (node == graph.root()) {
      return l;
    }

    node = to(ps.front());
    ps = graph.upEdges(node);
    l++;
  }
  return l;
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::inEdges(Edge edge) const {
  std::vector<Edge> prev = inEdges(from(edge));
  eraseValue(prev, flip(edge));
  return prev;
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::outEdges(Edge edge) const {
  std::vector<Edge> post = outEdges(to(edge));
  eraseValue(post, flip(edge));
  return post;
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::upEdges(Node node) const {
  std::map<Edge, Attribute> filtered_edges;
  std::copy_if(edges_.begin(), edges_.end(),
               std::inserter(filtered_edges, filtered_edges.end()),
               [node](const std::pair<Edge, Attribute> &p) {
                 return (from(p.first) == node) && isUpEdge(p.first);
               });

  return mapToVector(filtered_edges);
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::upEdges(Edge edge) const {
  return upEdges(to(edge));
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::downEdges(Node node) const {
  std::map<Edge, Attribute> filtered_edges;
  std::copy_if(edges_.begin(), edges_.end(),
               std::inserter(filtered_edges, filtered_edges.end()),
               [node](const std::pair<Edge, Attribute> &p) {
                 return (from(p.first) == node) && isDownEdge(p.first);
               });

  return mapToVector(filtered_edges);
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::downEdges(Edge edge) const {
  return downEdges(to(edge));
}

template <class Attribute>
std::vector<Edge> Graph<Attribute>::sortedEdges() const {
  std::map<Edge, Attribute> filtered_edges;
  /// copy up-edges, i.e. where from(edge) < to(edge)
  std::copy_if(
      edges_.begin(), edges_.end(),
      std::inserter(filtered_edges, filtered_edges.end()),
      [](const std::pair<Edge, Attribute> &p) { return (isUpEdge(p.first)); });
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

template <class Attribute>
index_t Graph<Attribute>::outIndex(Edge edge) const {
  return inIndex(flip(edge));
/*  Node node = from(edge);
  auto vec = outEdges(node);
  auto it = std::find(vec.begin(), vec.end(), edge);
  if (it != vec.end()) {
    auto index = std::distance(vec.begin(), it);
    return (index_t)index;
  } else {
    throw std::runtime_error("outindex: cannot find edge.");
    return 0;
  }*/
}

template <class Attribute>
index_t Graph<Attribute>::inIndex(Edge edge) const {
  Node node = to(edge);
  auto vec = inEdges(node);
  auto it = std::find(vec.begin(), vec.end(), edge);
  if (it != vec.end()) {
    auto index = std::distance(vec.begin(), it);
    return (index_t)index;
  } else {
    throw std::runtime_error("inindex: cannot find edge.");
    return 0;
  }
}

template <class Attribute> std::vector<Node> Graph<Attribute>::leaves() const {
  std::vector<Node> ls;
  for (const auto& edge : sortedEdges()) {
    const Node node = from(edge);
    if (node < 0) { ls.push_back(node);}
  }
  return ls;
}

template class Graph<tensorlib::Tensor>;
template class Graph<std::vector<long long>>;
template class Graph<std::any>;

} // namespace qutree

std::ostream &operator<<(std::ostream &os, const qutree::Edge &edge) {
  os << "{" << qutree::from(edge) << ", " << qutree::to(edge) << "}";
  return os;
}

std::ostream &operator<<(std::ostream &os, const qutree::Graph<> &graph) {
  os << "Nodes:\n";
  for (auto p : graph.nodes_) {
    std::any anything = graph.nodes_.at(p.first);
    os << p.first;
    if (anything.type() == typeid(std::string)) {
      os << "\t" << std::any_cast<std::string>(anything);
    }
    os << std::endl;
  }

  os << "Edges:\n";
  for (qutree::Edge e : graph.sortedEdges()) {
    std::any anything = graph.edges_.at(e);
    os << e;
    if (anything.type() == typeid(std::string)) {
      os << "\t" << std::any_cast<std::string>(anything);
    }

    os << std::endl;
  }

  return os;
}
