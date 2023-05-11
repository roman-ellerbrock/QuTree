#pragma once
#include "torch/torch.h"
#include "torch_util.h"

namespace qutree {

namespace tensorlib = torch;
using Tensor = torch::Tensor;
using index_t = long long;

namespace torch {

at::TensorOptions options();
} // namespace torch

} // namespace qutree