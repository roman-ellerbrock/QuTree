#include "GraphFactory.h"
#include "OperatorFactory.h"
#include "TensorNetwork.h"
#include "arithmetics.h"
#include "gtest/gtest.h"

using namespace qutree;
using namespace std;

TEST(Arithmetics, dotProduct) {
  Graph graph = balancedBinaryTree(4);
  NetworkShape shape = standardShape(graph, 3, 5);
  TensorNetwork tn = createTN(shape, CTN);
  TensorNetwork mt = createTN(shape, MN);
  dotProduct(mt, tn, tn);
}

TEST(Arithmetics, contract) {
  Graph graph = balancedBinaryTree(4);
  NetworkShape shape = standardShape(graph, 3, 5);
  TensorNetwork tn = createTN(shape, CTN);
  TensorNetwork mt = createTN(shape, MN);
  dotProduct(mt, tn, tn);
  ProductOperator P = identity(shape);
  TensorNetwork mt2 = matrixNetwork(shape, P);
  contract(mt2, tn, tn);
  // todo(Roman): only equal if qr implemented, so do it!
  for (auto e : shape.sortedEdges()) {
    //    ASSERT_NEAR(0., (mt.edges_[e]-mt2.edges_[e]).norm().item<double>(),
    //    1e-10);
  }
}

/*TEST(Arithmetics, qr) {
  Graph graph = balancedBinaryTree(4);
  NetworkShape shape = standardShape(graph, 3, 5);
  TensorNetwork tn = createTN(shape, CTN);
  TensorNetwork mt = createTN(shape, MN);
  dotProduct(mt, tn, tn);

  TensorNetwork umt = mt;
  TensorNetwork utn = qr(tn);
  dotProduct(umt, utn, utn);
  std::cout << umt << std::endl;
  std::cout << mt << std::endl;

  // even lowest mat is multiplied by underlying R
  Edge edge({4, 0});
  Tensor rho = mt.edges_[edge];
  Tensor urho = umt.edges_[edge];
  Tensor A = tn.nodes_[(to(edge))];
  std::cout << "A shape: " << A.sizes() << std::endl;
  auto qr = tensorlib::qr(A, 1);
  auto R = get<1>(qr);
  std::cout << edge << std::endl;
  std::cout << rho << std::endl;
  std::cout << urho << std::endl;
  rho = torch::mm(R.t(), rho);
  rho = torch::mm(rho, R);
  std::cout << rho << std::endl;
  urho = torch::mm(R.t(), urho);
  urho = torch::mm(urho, R);
  std::cout << urho << std::endl;
}*/

TEST(Arithmetics, two) {
  using namespace tensorlib;
  torch::Tensor A = torch::rand({2, 4, 3}, torch::options());
  torch::Tensor B = torch::rand({2, 3, 4}, torch::options());

  // manually calc original dot
  auto x = tensorlib::contractTensorTensor(A, A, 2);
  auto Bt = contractMatrixTensor(x, B, 1);
  x = contractTensorTensor(B, Bt, 2);
  //  std::cout << x << std::endl;

  // manually qr
  auto QR = tensorlib::qr(A, 2);
  auto Q = std::get<0>(QR);
  auto R = std::get<1>(QR);
  auto uA = Q;
  auto uB = contractMatrixTensor(R, B, 1);

  // recalc dot
  x = contractTensorTensor(uA, uA, 2);
  auto uBt = contractMatrixTensor(x, uB, 1);
  x = contractTensorTensor(uB, uB, 2);
  //  std::cout << x << std::endl;
}

TEST(Arithmetics, dot) {
  NetworkShape shape = testTree();
  TensorNetwork tn = createTN(shape, CTN);
  TensorNetwork mt = createTN(shape, MN);
  dotProduct(mt, tn, tn);

  // calculate manual
  using namespace tensorlib;
  auto TN = tn;
  auto MT = mt;
  Tensor A0 = TN.nodes_[0];
  Tensor A1 = TN.nodes_[1];
  Tensor A2 = TN.nodes_[2];
  MT.edges_[{0, 2}] = contractTensorTensor(A0, A0, 2);
  MT.edges_[{1, 2}] = contractTensorTensor(A1, A1, 1);

  auto xA2 = contractMatrixTensor(MT.edges_[{0, 2}], A2, 0);
  MT.edges_[{2, 1}] = contractTensorTensor(A2, xA2, 1);

  xA2 = contractMatrixTensor(MT.edges_[{1, 2}], A2, 1);
  MT.edges_[{2, 0}] = contractTensorTensor(A2, xA2, 0);

  for (Edge edge : MT.sortedEdges()) {
    auto M = MT.edges_[edge];
    auto m = mt.edges_[edge];
    ASSERT_NEAR(0., (m - M).norm().item<double>(), 1e-12);
  }
}

TEST(Arithmetics, qr) {
  NetworkShape shape = testTree();
  TensorNetwork tn = createTN(shape, CTN);

  auto utn = qr(tn);
  TensorNetwork umt = createTN(shape, MN);
  dotProduct(umt, utn, utn);

  // calculate manual
  using namespace tensorlib;
  auto TN = tn;
  auto MT = umt;
  Tensor &A2 = TN.nodes_[2];
  // qr on 0
  Tensor &A0 = TN.nodes_[0];
  auto QR = tensorlib::qr(A0, 2);
  auto Q = get<0>(QR);
  auto R = get<1>(QR);
  A0 = Q;
  A2 = contractMatrixTensor(R, A2, 0);

  // qr on 1
  Tensor &A1 = TN.nodes_[1];
  QR = tensorlib::qr(A1, 1);
  Q = get<0>(QR);
  R = get<1>(QR);
  A1 = Q;
  A2 = contractMatrixTensor(R, A2, 1);
  std::cout << TN << std::endl;
  std::cout << tn << std::endl;

  dotProduct(MT, TN, TN);

  for (Edge edge : MT.sortedEdges()) {
    auto M = MT.edges_[edge];
    auto m = umt.edges_[edge];
    ASSERT_NEAR(0., (m - M).norm().item<double>(), 1e-12);
  }
}

TEST(Arithmetics, qr2) {
  using namespace tensorlib;
  /// compare dot after qr with the one before
  NetworkShape shape = testTree();
  TensorNetwork tn = createTN(shape, CTN);

  // dot after qr
  auto utn = qr(tn);
  TensorNetwork umt = createTN(shape, MN);
  dotProduct(umt, utn, utn);

  // there
  auto TN = tn;
  auto MT = umt;
  // qr on 1
  Tensor &A1 = TN.nodes_[1];
  Tensor &A2 = TN.nodes_[2];
  auto QR = tensorlib::qr(A1, 1);
  auto Q = get<0>(QR);
  auto R = get<1>(QR);
  A1 = Q;
  A2 = contractMatrixTensor(R, A2, 1);

  dotProduct(MT, TN, TN);

  auto M = MT.edges_[{2, 1}];
  auto m = umt.edges_[{2, 1}];
  ASSERT_NEAR(0., (m - M).norm().item<double>(), 1e-12);
}