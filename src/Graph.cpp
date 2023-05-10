#include "Graph.h"

namespace qutree
{
    Node from(Edge e) { return e.first; }
    Node to(Edge e) { return e.second; }
    Edge flip(Edge e) { return {to(e), from(e)}; }

    std::ostream& operator<<(std::ostream& os, const Edge& e) {
        os << "{" << e.first << ", " << e.second << "}";
        return os;
    }

    template <class Attribute>
    std::vector<Edge> mapToVector(const std::map<Edge, Attribute> &m)
    {
        std::vector<Edge> v;
        std::transform(m.begin(), m.end(), std::back_inserter(v),
                       [](const std::pair<Edge, Attribute> &p)
                       { return p.first; });
        return v;
    }

    void eraseValue(std::vector<Edge> &edges, Edge e)
    {
        while (true)
        {
            auto it = std::find(edges.begin(), edges.end(), e);
            if (it == edges.end())
            {
                return;
            }
            edges.erase(it);
        }
    }

    bool isUpEdge(Edge e)
    {
        return (to(e) > from(e));
    }

    bool isDownEdge(Edge e)
    {
        return !isUpEdge(e);
    }
}


std::ostream& operator<<(std::ostream& os, const qutree::Edge& edge)
{
    os << "{" << qutree::from(edge) << ", " << qutree::to(edge) << "}";
    return os;
}

std::ostream& operator<<(std::ostream& os, const qutree::Graph<>& graph)
{
    os << "Leaves:\n";
    for (auto p : graph.leaves_)
    {
        qutree::Leaf l = p.first;
        os << l << std::endl;
    }

    os << "Nodes:\n";
    for (auto p : graph.nodes_)
    {
        os << p.first << std::endl;
    }
    os << "Edges:\n";
    for (qutree::Edge e : graph.sortedEdges())
    {
        os << e << std::endl;
    }
    return os;   
}
