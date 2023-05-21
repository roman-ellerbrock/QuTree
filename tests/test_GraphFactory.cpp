#include "GraphFactory.h"
#include "graph_util.h"
#include "gtest/gtest.h"
#include <ranges>

using namespace qutree;

TEST(GraphFactory, subgraph) {
  /**
   *       6
   *     /   \
   *   4       5
   *  / \     / \
   * 0   1   2   3
   * |   |   |   |
   *-1  -2  -3  -4
   *     X   X
   **/
  Graph graph = binary_4_graph();
  std::vector<Leaf> leaves({{-2, 1}, {-3, 2}});
  Graph sub = subgraph(graph, leaves);

  /// Nodes
  ASSERT_TRUE(sub.containsNode(1));
  ASSERT_TRUE(sub.containsNode(2));
  ASSERT_TRUE(sub.containsNode(4));
  ASSERT_TRUE(sub.containsNode(5));
  ASSERT_TRUE(sub.containsNode(6));

  ASSERT_FALSE(sub.containsNode(0));
  ASSERT_FALSE(sub.containsNode(3));

  /// Leaves
  ASSERT_TRUE(sub.containsEdge({-2, 1}));
  ASSERT_TRUE(sub.containsEdge({-3, 2}));

  ASSERT_FALSE(sub.containsEdge({-1, 0}));
  ASSERT_FALSE(sub.containsEdge({-4, 3}));

  /// Edges
  ASSERT_TRUE(sub.containsEdge({1, 4}));
  ASSERT_TRUE(sub.containsEdge({4, 6}));
  ASSERT_TRUE(sub.containsEdge({2, 5}));
  ASSERT_TRUE(sub.containsEdge({5, 6}));

  ASSERT_FALSE(sub.containsEdge({0, 4}));
  ASSERT_FALSE(sub.containsEdge({3, 5}));
}

TEST(GraphFactory, balancedBinaryTree) {
  /**
   *        4
   *      /   \
   *     3     2
   *   /   \
   *  0     1
   */
  Graph<> graph = balancedBinaryTree(3);
  ASSERT_TRUE(graph.containsEdge({-1, 0}));
  ASSERT_TRUE(graph.containsEdge({-2, 1}));
  ASSERT_TRUE(graph.containsEdge({-3, 2}));

  ASSERT_TRUE(graph.containsNode(3));
  ASSERT_TRUE(graph.containsNode(4));

  ASSERT_TRUE(graph.containsEdge({0, 3}));
  ASSERT_TRUE(graph.containsEdge({1, 3}));
  ASSERT_TRUE(graph.containsEdge({2, 4}));
  ASSERT_TRUE(graph.containsEdge({3, 4}));
  ASSERT_TRUE(graph.containsEdge(flip({0, 3})));
  ASSERT_TRUE(graph.containsEdge(flip({1, 3})));
  ASSERT_TRUE(graph.containsEdge(flip({2, 4})));
  ASSERT_TRUE(graph.containsEdge(flip({3, 4})));
}
