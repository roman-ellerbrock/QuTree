#include "gtest/gtest.h"
#include "TensorNetwork.h"
#include <ranges>
#include "graph_util.h"

using namespace qutree;


TEST (TensorNetwork, construct)
{
    Graph graph = binary_4_graph();
    NetworkShape shapes = createNetworkShape(graph, 4, 10);
    std::cout << shapes << std::endl;
}