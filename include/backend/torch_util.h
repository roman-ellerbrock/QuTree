#pragma once
#include "torch/torch.h"

namespace torch {

at::TensorOptions options();

// Wrapper lambda function for torch::eye
Tensor eyeWrapper(IntArrayRef size, TensorOptions options);

// Reshape a tensor to third-order tensor by selecting the k-th index
Tensor reshape3(Tensor A, long long k);

// contract a matrix with a 3rd order tensor
Tensor contract2_3(Tensor mat, Tensor A);

// contract two 3rd order tensors
Tensor contract3_3(Tensor A, Tensor B);

// contract matrix with tensor
Tensor contractMatrixTensor(Tensor mat, Tensor A, long long idx);

// contract tensor with tensor
Tensor contractTensorTensor(Tensor A, Tensor B, long long idx);

// qr decomposition for tensor
/**
 * \brief perform qr decomposition matrixTensor(R.T, Q, idx) = A
*/
std::tuple<Tensor, Tensor> qr(Tensor A, long long idx);

} // namespace torch
