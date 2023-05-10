#include "TensorNetwork.h"

namespace qutree {

NetworkShape createNetworkShape(const Graph<> &graph,
    int64_t bondDimension,
    int64_t leafDimension) 
{
    NetworkShape network_shape(graph);

    for (auto& leafpair : network_shape.leaves_)
    {
        torch::IntArrayRef& ref = leafpair.second;
        ref = std::vector<int64_t>({leafDimension});
    }

    for (auto& edgepair : network_shape.edges_)
    {
        torch::IntArrayRef& ref = edgepair.second;
        ref = std::vector<int64_t>({bondDimension, bondDimension});
    }

    for (auto& nodepair : network_shape.nodes_)
    {
        torch::IntArrayRef& ref = nodepair.second;
        Node node = nodepair.first;
        auto nNeighbors = graph.outEdges(node).size();
        ref = std::vector<int64_t>(nNeighbors, bondDimension);
    }

    return network_shape;
}

TensorNetwork createTensorNetwork(const Graph<> &graph) {

  TensorNetwork psi(graph);

  return psi;
}

} // namespace qutree