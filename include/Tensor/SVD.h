//
// Created by Roman Ellerbrock on 11/18/21.
//

#ifndef SVD_H
#define SVD_H
#include "Tensor.h"

template<typename T>
class SVD: public tuple<Tensor<T>, Tensor<T>, Tensord> {
public:
	SVD() = default;
	~SVD() = default;

	explicit SVD(const TensorShape& shape);

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


#endif //SVD_H
