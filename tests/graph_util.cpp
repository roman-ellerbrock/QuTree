#include "graph_util.h"

namespace qutree {

std::map<Node, std::any> create_nodes() {
  return {{0, ""}, {1, ""}, {2, "A"}, {3, ""}, {4, ""}, {5, "B"}, {6, "root"}};
}

std::map<Edge, std::any> create_edges() {
  return {{{0, 2}, ""}, {{1, 2}, ""}, {{2, 6}, "A-r"},
          {{3, 5}, ""}, {{4, 5}, ""}, {{5, 6}, "B-r"},

          {{2, 0}, ""}, {{2, 1}, ""}, {{6, 2}, "A-r"},
          {{5, 3}, ""}, {{5, 4}, ""}, {{6, 5}, "B-r"}};
}

std::map<Leaf, std::any> create_leaves() {
  return {{{0, 0}, ""}, {{1, 1}, ""}, {{2, 3}, ""}, {{3, 4}, ""}};
}

Graph<> binary_4_graph() {
  return Graph(create_nodes(), create_edges(), create_leaves());
}

} // namespace qutree
