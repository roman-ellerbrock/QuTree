//
// Created by Roman Ellerbrock on 11/18/21.
//

#ifndef SVD_H
#define SVD_H
#include "Tensor.h"
#include <functional>

template<typename T>
class SVD: public tuple<Tensor<T>, Tensor<T>, Tensord> {
public:
	SVD() = default;
	~SVD() = default;

	explicit SVD(const TensorShape& shape)
		: SVD(shape, shape.lastIdx()) {
	}

	SVD(const TensorShape& shape, size_t k);

	explicit SVD(const tuple<Tensor<T>, Tensor<T>, Tensord>& tp) {
		U() = get<0>(tp);
		VT() = get<1>(tp);
		sigma() = get<2>(tp);
	}

	[[nodiscard]] const Tensor<T>& U() const { return get<0>(*this); }

	[[nodiscard]] const Tensor<T>& VT() const { return get<1>(*this); }

	[[nodiscard]] const Tensor<double>& sigma() const { return get<2>(*this); }

	Tensor<T>& U() { return get<0>(*this); }

	Tensor<T>& VT() { return get<1>(*this); }

	Tensor<double>& sigma() { return get<2>(*this); }
};

typedef SVD<complex<double>> SVDcd;

typedef SVD<double> SVDd;

template<typename T>
class SpectralDecomposition: public pair<Tensor<T>, Tensord> {
public:
	SpectralDecomposition() = default;
	~SpectralDecomposition() = default;

	explicit SpectralDecomposition(const tuple<Tensor<T>, Tensord>& tp) {
		U() = get<0>(tp);
		ev() = get<1>(tp);
	}

	SpectralDecomposition(const TensorShape& shape) {
		assert(shape.order() == 2);
		assert(shape[0] == shape[1]);
		Tensor<T> U({shape[0], shape[1]});
		Tensord ev({shape[0]});
		*this = SpectralDecomposition<T>({U, ev});
	}

	Tensor<T>& U() { return get<0>(*this); }

	Tensor<double>& ev() { return get<1>(*this); }

	[[nodiscard]] const Tensor<T>& U() const { return get<0>(*this); }

	[[nodiscard]] const Tensor<double>& ev() const { return get<1>(*this); }
};

typedef SpectralDecomposition<complex<double>> SpectralDecompositioncd;

typedef SpectralDecomposition<double> SpectralDecompositiond;

template<typename T>
class FactorizedGe: public tuple<Tensor<T>, Tensor<T>, Tensor<T>> {
public:
	FactorizedGe() = default;
	~FactorizedGe() = default;

	explicit FactorizedGe(const TensorShape& shape) {
		assert(shape.order() == 2);
		size_t m = shape[0];
		size_t n = shape[1];
		assert(n == m);
		VR() = Tensor<T>({m, n});
		VL() = Tensor<T>({m, n});
		ev() = Tensor<T>({n});

	}

	explicit FactorizedGe(const tuple<Tensor<T>, Tensor<T>, Tensor<T>>& tp) {
		VR() = get<0>(tp);
		VL() = get<1>(tp);
		ev() = get<2>(tp);
	}

	Tensor<T>& VL() { return get<0>(*this); }

	[[nodiscard]] const Tensor<T>& VL() const { return get<0>(*this); }

	Tensor<T>& ev() { return get<1>(*this); }

	[[nodiscard]] const Tensor<T>& ev() const { return get<1>(*this); }

	Tensor<T>& VR() { return get<2>(*this); }

	[[nodiscard]] const Tensor<T>& VR() const { return get<2>(*this); }

};

#endif //SVD_H
