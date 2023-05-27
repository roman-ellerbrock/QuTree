#pragma once
#include "torch/torch.h"
#include "torch_util.h"
#include <complex>

namespace qutree {

namespace tensorlib = torch;

using Tensor = torch::Tensor;
using index_t = long long;

using real_t = double;
using complex_t = std::complex<double>;

using scalar_t = complex_t;

} // namespace qutree

namespace torch {

at::TensorOptions options();
at::TensorOptions real_options();
at::TensorOptions complex_options();
} // namespace torch
