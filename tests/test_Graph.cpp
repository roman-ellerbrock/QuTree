#include "gtest/gtest.h"
#include "Graph.h"
#include <ranges>
#include "graph_util.h"

using namespace qutree;


TEST(Graph, initialize)
{
    Graph graph = binary_4_graph();
    EXPECT_EQ(7, graph.nodes_.size());
    EXPECT_EQ(12, graph.edges_.size());
    EXPECT_EQ(4, graph.leaves_.size());
}

TEST(Graph, inEdges)
{
    Graph graph = binary_4_graph();
    std::vector<Edge> es = graph.inEdges(2);
    EXPECT_EQ(Edge({0, 2}), es[0]);
    EXPECT_EQ(Edge({1, 2}), es[1]);
    EXPECT_EQ(Edge({6, 2}), es[2]);
}

TEST(Graph, outEdges)
{
    Graph graph = binary_4_graph();
    std::vector<Edge> es = graph.outEdges(2);
    EXPECT_EQ(Edge({2, 0}), es[0]);
    EXPECT_EQ(Edge({2, 1}), es[1]);
    EXPECT_EQ(Edge({2, 6}), es[2]);
}

TEST(Graph, inEdges_edge)
{
    Graph graph = binary_4_graph();
    std::vector<Edge> es = graph.inEdges({2, 0});
    EXPECT_EQ(2, es.size());
    EXPECT_EQ(Edge({1, 2}), es[0]);
    EXPECT_EQ(Edge({6, 2}), es[1]);

    es = graph.inEdges({2, 6});
    EXPECT_EQ(2, es.size());
    EXPECT_EQ(Edge({0, 2}), es[0]);
    EXPECT_EQ(Edge({1, 2}), es[1]);
}

TEST(Graph, outEdges_edge)
{
    Graph graph = binary_4_graph();
    std::vector<Edge> es = graph.outEdges({0, 2});
    EXPECT_EQ(2, es.size());
    EXPECT_EQ(Edge({2, 1}), es[0]);
    EXPECT_EQ(Edge({2, 6}), es[1]);

    es = graph.outEdges({6, 2});
    EXPECT_EQ(2, es.size());
    EXPECT_EQ(Edge({2, 0}), es[0]);
    EXPECT_EQ(Edge({2, 1}), es[1]);
}

TEST(Graph, empty)
{
    Graph graph;
    ASSERT_TRUE(graph.empty());
    graph = binary_4_graph();
    ASSERT_FALSE(graph.empty());
}

TEST(Graph, clear)
{
    Graph graph = binary_4_graph();
    ASSERT_FALSE(graph.empty());
    graph.clear();
    ASSERT_TRUE(graph.empty());
}

TEST(Graph, containsNode)
{
    Graph graph = binary_4_graph();
    ASSERT_TRUE(graph.containsNode(0));
    ASSERT_FALSE(graph.containsNode(7));
}

TEST(Graph, containsEdge)
{
    Graph graph = binary_4_graph();
    ASSERT_TRUE(graph.containsEdge({0, 2}));
    ASSERT_FALSE(graph.containsEdge({0, 6}));
}

TEST(Graph, containsLeaf)
{
    Graph graph = binary_4_graph();
    ASSERT_TRUE(graph.containsLeaf({0, 0}));
    ASSERT_FALSE(graph.containsLeaf({0, 1}));
}

TEST(Graph, upEdges)
{
    Graph graph = binary_4_graph();
    auto es = graph.upEdges({0, 2});
    EXPECT_EQ(1, es.size());
    EXPECT_EQ(Edge({2, 6}), es[0]);
}

TEST(Graph, downEdges)
{
    Graph graph = binary_4_graph();
    auto es = graph.downEdges({6, 2});
    EXPECT_EQ(2, es.size());
    EXPECT_EQ(Edge({2, 0}), es[0]);
    EXPECT_EQ(Edge({2, 1}), es[1]);
}

TEST(Graph, iterateNodes)
{
    Graph graph = binary_4_graph();
    std::vector<Node> nodes = {0, 1, 2, 3, 4, 5, 6};
    size_t idx = 0;
    for (auto const &[node, name] : graph.nodes_)
    {
        EXPECT_EQ(nodes[idx++], node);
    }

    idx = nodes.size() - 1;
    for (auto it = graph.nodes_.rbegin(); it != graph.nodes_.rend(); ++it)
    {
        const auto &kv = *it;
        Node node = kv.first;
        EXPECT_EQ(nodes[idx--], node);
    }
}

TEST(Graph, iterateEdges)
{
    Graph graph = binary_4_graph();
    std::vector<Edge> es(
        {{0, 2},
         {1, 2},
         {2, 6},
         {3, 5},
         {4, 5},
         {5, 6},
         {6, 5},
         {6, 2},
         {5, 4},
         {5, 3},
         {2, 1},
         {2, 0}
         });
    size_t idx = 0;
    auto sorted_edges = graph.sortedEdges();
    EXPECT_EQ(sorted_edges.size(), es.size());
    for (const auto &edge : sorted_edges)
    {
        EXPECT_EQ(es[idx++], edge);
    }
}

TEST(Graph, layer)
{
    Graph graph = binary_4_graph();
    std::vector<Node> nodes = {0, 1, 2, 3, 4, 5, 6};
    ASSERT_EQ(2, qutree::layer(0, graph));
    ASSERT_EQ(2, qutree::layer(1, graph));
    ASSERT_EQ(2, qutree::layer(3, graph));
    ASSERT_EQ(2, qutree::layer(4, graph));

    ASSERT_EQ(1, qutree::layer(2, graph));
    ASSERT_EQ(1, qutree::layer(5, graph));

    ASSERT_EQ(0, qutree::layer(6, graph));
}
