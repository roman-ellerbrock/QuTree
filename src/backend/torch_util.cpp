#include "backend/Tensor.h"

namespace torch {

at::TensorOptions options() {
  return at::TensorOptions().dtype(torch::kDouble);
}

at::TensorOptions real_options() {
  return at::TensorOptions().dtype(torch::kDouble);
}

at::TensorOptions complex_options() {
  return at::TensorOptions().dtype(torch::kComplexDouble);
}

// Wrapper lambda function for torch::eye
Tensor eyeWrapper(IntArrayRef size, TensorOptions options) {
  auto n = size[0]; // Assuming a square matrix, so size[0] = size[1]
  return torch::eye(n, options);
};

Tensor reshape3(Tensor A, long long k) {
  auto shape = A.sizes();
  auto a = 1ll;
  for (auto i = 0ll; i < k; ++i) {
    a *= shape[i];
  }

  auto b = shape[k];

  auto c = 1ll;
  for (auto i = (k + 1); i < shape.size(); ++i) {
    c *= shape[i];
  }
  return A.reshape({a, b, c});
}

Tensor contract2_3(Tensor mat, Tensor A) {
  auto shape = A.sizes();
  auto a = shape[0];
  auto b = shape[1];
  auto c = shape[2];

  A = A.transpose(0, 1);
  A = A.reshape({b, a * c});
  Tensor B = torch::matmul(mat, A);
  B = B.reshape({b, a, c});
  B = B.transpose(0, 1);
  return B;
}

Tensor contract3_3(Tensor A, Tensor B) {
  /// prepare A
  auto shape = A.sizes();
  auto a = shape[0];
  auto b = shape[1];
  auto c = shape[2];
  A = A.transpose(0, 1);
  A = A.reshape({b, a * c});

  /// prepare B
  shape = B.sizes();
  a = shape[0];
  b = shape[1];
  c = shape[2];
  B = B.transpose(1, 2);
  B = B.reshape({a * c, b});

  return torch::mm(A.conj(), B);
}

Tensor contractMatrixTensor(Tensor mat, Tensor A, long long idx) {
  try {
    Tensor B;
    auto shape = A.sizes();
    A = torch::reshape3(A, idx);
    B = torch::contract2_3(mat, A);
    B = B.reshape(shape);
    return B;
  } catch (const std::exception &e) {
    std::cerr << "std::exception in matrix-tensor contraction:" << std::endl;
    std::cerr << e.what() << std::endl;
    return torch::Tensor();
  } catch (...) {
    std::cerr << "Unknown exception in matrix-tensor contraction." << std::endl;
    return torch::Tensor();
  }
}

Tensor contractTensorTensor(Tensor A, Tensor B, long long idx) {
  try {
    auto shapeA = A.sizes();
    A = torch::reshape3(A, idx);
    auto shapeB = B.sizes();
    B = torch::reshape3(B, idx);
    auto c = contract3_3(A, B);
    return c;
  } catch (const std::exception &e) {
    std::cerr << "std::exception in tensor-tensor contraction:" << std::endl;
    std::cerr << e.what() << std::endl;
    return torch::Tensor();
  } catch (...) {
    std::cerr << "Unknown exception in tensor-tensor contraction." << std::endl;
    return torch::Tensor();
  }
}

std::tuple<Tensor, Tensor> qr(Tensor A, long long idx) {
  // Q: tensor -> matrix
  auto shape = A.sizes();
  // -> {a, c, b}
  A = torch::reshape3(A, idx);
  A = A.transpose(1, 2);
  std::vector<long long> shape3(A.sizes().begin(), A.sizes().end());
  // -> {a * c, b}
  A = A.reshape({shape3[0] * shape3[1], shape3[2]});
  if (shape3[0] * shape3[1] < shape3[2]) {
    throw std::runtime_error("QR dimension violation: wrong rectangular shape.");
  }

  // do qr
  auto qr = at::linalg_qr(A, "reduced");

  // Q: matrix -> Tensor
  // -> {a, c, b}
  auto& Q = get<0>(qr);
  Q = Q.reshape(shape3);
  // -> {a, b, c}
  Q = Q.transpose(1, 2);
  // -> tensor
  Q = Q.reshape(shape);
  return qr;
}

} // namespace torch
