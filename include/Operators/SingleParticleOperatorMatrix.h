//
// Created by Roman Ellerbrock on 2020-01-24.
//

#ifndef SINGLEPARTICLEOPERATORMATRIX_H
#define SINGLEPARTICLEOPERATORMATRIX_H
#include "SingleParticleOperator.h"
#include "FactorMatrix.h"

template<typename T>
class SingleParticleOperatorMatrix: SPO<T> {
public:
	SingleParticleOperatorMatrix() = default;

	explicit SingleParticleOperatorMatrix(FactorMatrix<T> h);

	~SingleParticleOperatorMatrix() = default;

	virtual void Apply(const PrimitiveBasis& grid, Tensor<T>& hAcoeff,
		const Tensor<T>& Acoeff) const;
private:
	FactorMatrix<T> h_;
};

template<typename T>
using SPOM = SingleParticleOperatorMatrix<T>;

typedef SPOM<complex<double>> SPOMcd;

typedef SPOM<double> SPOMd;

#endif //SINGLEPARTICLEOPERATORMATRIX_H
