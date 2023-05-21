#include "graph_util.h"

namespace qutree {

std::map<Node, std::any> create_nodes() {
  return {{0, ""}, {1, ""}, {2, ""}, {3, ""}, {4, "A"}, {5, "B"}, {6, "root"}};
}

std::map<Edge, std::any> create_edges() {
  return {{{-1, 0}, ""},   {{-2, 1}, ""},   {{-3, 2}, ""}, {{-4, 3}, ""},

          {{0, 4}, ""},    {{1, 4}, ""},    {{2, 5}, ""},  {{3, 5}, ""},
          {{4, 6}, "A-r"}, {{5, 6}, "B-r"},

          {{4, 0}, ""},    {{4, 1}, ""},    {{5, 2}, ""},  {{5, 3}, ""},
          {{6, 4}, "A-r"}, {{6, 5}, "B-r"}};
}

Graph<> binary_4_graph() { return Graph(create_nodes(), create_edges()); }

} // namespace qutree
