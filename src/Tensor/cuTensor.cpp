#include "cuTensor.h"
#include "TensorBLAS1.h"
#include <lapack.hh>

using f =  float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;

//using namespace polymorphic;
//template class Tensor<float, cuMemory>;
//template class Tensor<double, cuMemory>;
//template class Tensor<complex<double>, cuMemory>;

using namespace polymorphic;

template<typename T, template <typename> class Mem, template <typename> class oMem>
Tensor<T, Mem> transfer(const Tensor<T, oMem>& src) {
	Tensor<T, Mem> dest(src.shape_);
	dest.mem() = src.mem();
	return dest;
}

template Tensor<d, hostMemory> transfer(const Tensor<d, cuMemory>& src);
template Tensor<f, hostMemory> transfer(const Tensor<f, cuMemory>& src);
template Tensor<cf, hostMemory> transfer(const Tensor<cf, cuMemory>& src);
template Tensor<cd, hostMemory> transfer(const Tensor<cd, cuMemory>& src);

template Tensor<f, cuMemory> transfer(const Tensor<f, hostMemory>& src);
template Tensor<d, cuMemory> transfer(const Tensor<d, hostMemory>& src);
template Tensor<cd, cuMemory> transfer(const Tensor<cd, hostMemory>& src);
template Tensor<cf, cuMemory> transfer(const Tensor<cf, hostMemory>& src);
