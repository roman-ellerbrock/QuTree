#include "OperatorFactory.h"
#include <gtest/gtest.h>
#include "GraphFactory.h"

using namespace qutree;
using namespace std;

TEST(OperatorFactory, pauliX) {
  auto x = pauli::X();
  ASSERT_NEAR(0., x[0][0].item<double>(), 1e-12);
  ASSERT_NEAR(1., x[1][0].item<double>(), 1e-12);
  ASSERT_NEAR(1., x[0][1].item<double>(), 1e-12);
  ASSERT_NEAR(0., x[1][1].item<double>(), 1e-12);
}

TEST(OperatorFactory, pauliY) {
  auto x = tensorlib::view_as_real(pauli::Y());
  ASSERT_NEAR(0., x[0][0][0].item<double>(), 1e-12);
  ASSERT_NEAR(0., x[1][0][0].item<double>(), 1e-12);
  ASSERT_NEAR(0., x[0][1][0].item<double>(), 1e-12);
  ASSERT_NEAR(0., x[1][1][0].item<double>(), 1e-12);

  ASSERT_NEAR(0., x[0][0][1].item<double>(), 1e-12);
  ASSERT_NEAR(1., x[1][0][1].item<double>(), 1e-12);
  ASSERT_NEAR(-1., x[0][1][1].item<double>(), 1e-12);
  ASSERT_NEAR(0., x[1][1][1].item<double>(), 1e-12);
}

TEST(OperatorFactory, pauliZ) {
  auto x = pauli::Z();
  ASSERT_NEAR(1., x[0][0].item<double>(), 1e-12);
  ASSERT_NEAR(0., x[1][0].item<double>(), 1e-12);
  ASSERT_NEAR(0., x[0][1].item<double>(), 1e-12);
  ASSERT_NEAR(-1., x[1][1].item<double>(), 1e-12);
}

TEST(OperatorFactory, Ising) {
  auto H = transversalFieldIsing(1., 2., 10);
  ASSERT_EQ(H.size(), 20);
  ASSERT_NEAR(real(H[0].first), -1., 1e-12);
  auto P = H[0].second;
  ASSERT_EQ(2, P.size());
  ASSERT_NEAR(real(H[1].first), -2., 1e-12);
  P = H[1].second;
  ASSERT_EQ(1, P.size());
}

TEST(OperatorFactory, Identity) {
  Graph graph = balancedBinaryTree(4);
  NetworkShape shape = standardShape(graph, 3, 5);
  ProductOperator P = identity(shape);
  ASSERT_EQ(P.size(), 4);
  ASSERT_NEAR(0., (P[-1] - torch::eye(5)).norm().item<double>(), 1e-12);
  ASSERT_NEAR(0., (P[-2] - torch::eye(5)).norm().item<double>(), 1e-12);
  ASSERT_NEAR(0., (P[-3] - torch::eye(5)).norm().item<double>(), 1e-12);
  ASSERT_NEAR(0., (P[-4] - torch::eye(5)).norm().item<double>(), 1e-12);
}
