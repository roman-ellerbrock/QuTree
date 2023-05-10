#include "GraphFactory.h"
#include "graph_util.h"
#include "gtest/gtest.h"
#include <ranges>

using namespace qutree;

TEST(GraphFactory, subgraph) {
  Graph graph = binary_4_graph();
  std::vector<Leaf> leaves({{1, 1}, {2, 3}});
  Graph sub = subgraph(graph, leaves);

  /// Leaves
  ASSERT_TRUE(sub.containsLeaf({1, 1}));
  ASSERT_TRUE(sub.containsLeaf({2, 3}));

  ASSERT_TRUE(!sub.containsLeaf({0, 0}));
  ASSERT_TRUE(!sub.containsLeaf({3, 4}));

  /// Nodes
  ASSERT_TRUE(sub.containsNode(1));
  ASSERT_TRUE(sub.containsNode(2));
  ASSERT_TRUE(sub.containsNode(3));
  ASSERT_TRUE(sub.containsNode(5));
  ASSERT_TRUE(sub.containsNode(6));

  /// Edges
  ASSERT_TRUE(sub.containsEdge({1, 2}));
  ASSERT_TRUE(sub.containsEdge({2, 6}));
  ASSERT_TRUE(sub.containsEdge({3, 5}));
  ASSERT_TRUE(sub.containsEdge({5, 6}));

  /// edges that should not be contained
  ASSERT_TRUE(!sub.containsEdge({0, 2}));
  ASSERT_TRUE(!sub.containsEdge({4, 5}));
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
  ASSERT_TRUE(graph.containsLeaf({0, 0}));
  ASSERT_TRUE(graph.containsLeaf({1, 1}));
  ASSERT_TRUE(graph.containsLeaf({2, 2}));

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

