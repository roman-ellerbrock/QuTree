#include "OperatorFactory.h"
#include <complex>

namespace qutree {

namespace pauli {
Tensor X() {
  return tensorlib::tensor({{0., 1.}, {1., 0.}}).to(tensorlib::kFloat64);
}

Tensor Y() {
  auto real_part = torch::tensor({{0., 0.},{0., 0.}}, tensorlib::kFloat64);
  auto imag_part = torch::tensor({{0., -1.},{1., 0.}}, tensorlib::kFloat64);
  return torch::complex(real_part, imag_part);
}

Tensor Z() {
  return tensorlib::tensor({{1., 0.}, {0., -1.}}).to(tensorlib::kFloat64);
}
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