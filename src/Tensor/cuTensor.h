#pragma once
#include "Tensor/Tensor.h"
#include "cuMemory.h"

template <typename T>
using cuTensor = Tensor<T, polymorphic::cuMemory>;

using cuTensorf = cuTensor<float>;
using cuTensord = cuTensor<double>;
using cuTensorcd = cuTensor<complex<double>>;

template<typename T>
void gemm(cuTensor<T>& c, const cuTensor<T>& a, const cuTensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b,
	blas::Queue& queue);

template<typename T, template <typename> class Mem, template <typename> class oMem>
Tensor<T, Mem> transfer(const Tensor<T, oMem>& src);

constexpr auto transferToGPUf = transfer< float, polymorphic::cuMemory, polymorphic::hostMemory>;
constexpr auto transferToGPUd = transfer< double, polymorphic::cuMemory, polymorphic::hostMemory>;
constexpr auto transferToGPUcf = transfer< complex<float>, polymorphic::cuMemory, polymorphic::hostMemory>;
constexpr auto transferToGPUcd = transfer< complex<double>, polymorphic::cuMemory, polymorphic::hostMemory>;

constexpr auto transferFromGPUf = transfer< float, polymorphic::hostMemory, polymorphic::cuMemory>;
constexpr auto transferFromGPUd = transfer< double, polymorphic::hostMemory, polymorphic::cuMemory>;
constexpr auto transferFromGPUcf = transfer< complex<float>, polymorphic::hostMemory, polymorphic::cuMemory>;
constexpr auto transferFromGPUcd = transfer< complex<double>, polymorphic::hostMemory, polymorphic::cuMemory>;