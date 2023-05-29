#include "OperatorFactory.h"
#include <complex>

namespace qutree {

namespace pauli {
Tensor X() {
  return tensorlib::tensor({{0., 1.}, {1., 0.}}, tensorlib::kFloat64);
}

Tensor Y() {
  auto real_part = torch::tensor({{0., 0.},{0., 0.}}, tensorlib::kFloat64);
  auto imag_part = torch::tensor({{0., -1.},{1., 0.}}, tensorlib::kFloat64);
  return torch::complex(real_part, imag_part);
}

Tensor Z() {
  return tensorlib::tensor({{1., 0.}, {0., -1.}}, tensorlib::kFloat64);
}
}

Tensor identity(index_t n) {
  return tensorlib::eye(n, tensorlib::options());
}

ProductOperator identity(const NetworkShape& shapes) {
  ProductOperator P;
  std::vector<Node> leaves = shapes.leaves();
  for (Node leaf : leaves) {
    Edge leafEdge = shapes.leafEdge(leaf);
    auto shape = shapes.edges_.at(leafEdge);
    index_t dim = shape.front();
    P[leaf] = identity(dim);
  }
  return P;
}

SumOfProducts transversalFieldIsing(double J, double g, index_t N) {
  SumOfProducts H;

  for (index_t i = 0; i < N; ++i) {
    ProductOperator P;
    P[i] = pauli::Z();
    P[(i + 1) % N] = pauli::Z();
    H.push_back({-J, P});

    ProductOperator Q;
    Q[i] = pauli::X();
    H.push_back({-g, Q});
  }

  return H;
}

} // namespace qutree