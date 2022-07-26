//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Tensor/Tensor.hpp"
#include <tuple>

typedef double d;
typedef float f;
typedef complex<double> cd;
typedef complex<float> cf;



template class Tensor<f>;
template class Tensor<d>;
template class Tensor<cf>;
template class Tensor<cd>;

template class Tensor<int64_t>;

template Tensor<f> randomGen(const TensorShape& shape);
template Tensor<d> randomGen(const TensorShape& shape);
template Tensor<cf> randomGen(const TensorShape& shape);
template Tensor<cd> randomGen(const TensorShape& shape);

// https://stackoverflow.com/questions/50338955/how-can-i-concisely-write-a-lot-of-explicit-function-template-instantiations
template<typename... Ts>
auto instantiateTensorFactories() {
    static auto funcs = std::tuple_cat(std::make_tuple(
        arange<Ts>,
        random<Ts>,
        identity<Ts>,
        delta<Ts>
    )...);

    return &funcs;
}

template auto instantiateTensorFactories<f, d, cf, cd>();
