#include "GraphFactory.h"

namespace qutree
{
    bool contains(const std::vector<Edge> &edges, Edge edge)
    {
        for (Edge e : edges)
        {
            if (e == edge)
            {
                return true;
            }
        }
        return false;
    }

    template <class Attribute>
    void addToGraph(Graph<Attribute>& subGraph, const Graph<Attribute>& graph, Edge edge) {
        std::vector<Edge> edges = graph.upEdges(edge);
        for (Edge e : edges) {
            subGraph.edges_[e] = graph.edges_.at(e);
            subGraph.nodes_[to(e)] = graph.edges_.at(e);
            addToGraph(subGraph, graph, e);
        }
    }

    /**
     * \brief Select smallest subgraph that connects all leaves
     **/
    template <class Attribute>
    Graph<Attribute> subgraph(const Graph<Attribute> &graph, const std::vector<Leaf> &leaves)
    {
        Graph<Attribute> subgraph;
        for (Leaf leaf : leaves)
        {
            subgraph.leaves_[leaf] = graph.leaves_.at(leaf);
            subgraph.nodes_[to(leaf)] = graph.nodes_.at(to(leaf));
            addToGraph(subgraph, graph, leaf);
        }
        return subgraph;
    }

    template Graph<std::any> subgraph(const Graph<std::any> &graph, const std::vector<Leaf> &leaves);
}
