#include "backend/Tensor.h"

namespace torch {

at::TensorOptions options() { return at::TensorOptions().dtype(torch::kDouble); }

} // namespace torch