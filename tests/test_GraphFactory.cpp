#include "gtest/gtest.h"
#include "GraphFactory.h"
#include <ranges>
#include "graph_util.h"

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


