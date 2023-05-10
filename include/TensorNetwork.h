#include "Graph.h"
#include "backend/Tensor.h"

namespace qutree {

using TensorNetwork = Graph<Tensor>;
using NetworkShape = Graph<tensorlib::IntArrayRef>;

NetworkShape createNetworkShape(const Graph<> &graph, int64_t bondDimension,
                                int64_t leafDimension);
TensorNetwork createTensorNetwork(const Graph<> &graph);

} // namespace qutree
