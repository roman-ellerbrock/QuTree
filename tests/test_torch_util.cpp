#include "backend/torch_util.h"
#include "gtest/gtest.h"

TEST(torch_util, reshape3) {
  torch::Tensor A = torch::rand({2, 3, 4, 5, 6});
  auto B = torch::reshape3(A, 2);
  ASSERT_EQ(std::vector<long long>({6, 4, 30}), B.sizes());
}

TEST(torch_util, contract2_3) {
  torch::Tensor A = torch::arange(8, torch::options());
  A = A.reshape({2, 2, 2});
  torch::Tensor mat = torch::arange(4, torch::options());
  mat = mat.reshape({2, 2});

  auto B = torch::contract2_3(mat, A);
  ASSERT_EQ(std::vector<long long>({2, 2, 2}), B.sizes());

  double eps = 1e-12;
  B = B.reshape({8});
  ASSERT_NEAR(2., B[0].item<double>(), eps);
  ASSERT_NEAR(3., B[1].item<double>(), eps);
  ASSERT_NEAR(6., B[2].item<double>(), eps);
  ASSERT_NEAR(11., B[3].item<double>(), eps);

  ASSERT_NEAR(6., B[4].item<double>(), eps);
  ASSERT_NEAR(7., B[5].item<double>(), eps);
  ASSERT_NEAR(26., B[6].item<double>(), eps);
  ASSERT_NEAR(31., B[7].item<double>(), eps);
}

TEST(torch_util, contract3_3) {

  torch::Tensor A = torch::arange(8, torch::options());
  A = A.reshape({2, 2, 2});
  torch::Tensor B = A;

  auto mat = torch::contract3_3(A, B);
  ASSERT_EQ(std::vector<long long>({2, 2}), mat.sizes());

  double eps = 1e-12;
  ASSERT_NEAR(42., mat[0][0].item<double>(), eps);
  ASSERT_NEAR(62., mat[1][0].item<double>(), eps);
  ASSERT_NEAR(62., mat[0][1].item<double>(), eps);
  ASSERT_NEAR(98., mat[1][1].item<double>(), eps);
}

TEST(torch_util, matrixTensor) {
  torch::Tensor A = torch::arange(8, torch::options());
  A = A.reshape({2, 1, 2, 1, 2});
  torch::Tensor mat = torch::arange(4, torch::options());
  mat = mat.reshape({2, 2});

  long long k = 2;
  auto B = torch::contractMatrixTensor(mat, A, k);
  ASSERT_EQ(std::vector<long long>({2, 1, 2, 1, 2}), B.sizes());

  double eps = 1e-12;
  B = B.reshape({8});
  ASSERT_NEAR(2., B[0].item<double>(), eps);
  ASSERT_NEAR(3., B[1].item<double>(), eps);
  ASSERT_NEAR(6., B[2].item<double>(), eps);
  ASSERT_NEAR(11., B[3].item<double>(), eps);

  ASSERT_NEAR(6., B[4].item<double>(), eps);
  ASSERT_NEAR(7., B[5].item<double>(), eps);
  ASSERT_NEAR(26., B[6].item<double>(), eps);
  ASSERT_NEAR(31., B[7].item<double>(), eps);
}

TEST(torch_util, contractTensorTensor) {

  torch::Tensor A = torch::arange(8, torch::options());
  A = A.reshape({2, 1, 2, 1, 2});
  torch::Tensor B = A;

  auto mat = torch::contractTensorTensor(A, B, 2);
  ASSERT_EQ(std::vector<long long>({2, 2}), mat.sizes());

  double eps = 1e-12;
  ASSERT_NEAR(42., mat[0][0].item<double>(), eps);
  ASSERT_NEAR(62., mat[1][0].item<double>(), eps);
  ASSERT_NEAR(62., mat[0][1].item<double>(), eps);
  ASSERT_NEAR(98., mat[1][1].item<double>(), eps);
}
