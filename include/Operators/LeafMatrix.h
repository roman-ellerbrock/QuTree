//
// Created by Roman Ellerbrock on 2020-01-24.
//
#ifndef LEAFMATRIX_H
#define LEAFMATRIX_H
#include "LeafOperator.h"
#include "Core/FactorMatrix.h"

template<typename T>
class LeafMatrix: public LeafOperator<T>
	/**
	 * \class LeafMatrix
	 * \ingroup Operators
	 * \brief This class allows to create SPOs from (factor) Matrices.
	 */
{
public:
	LeafMatrix() = default;

	explicit LeafMatrix(FactorMatrix<T> h);

	explicit LeafMatrix(Matrix<T> h);

	~LeafMatrix() = default;

	virtual void Apply(const PrimitiveBasis& grid, Tensor<T>& hAcoeff,
		const Tensor<T>& Acoeff) const override;

	size_t& Mode() { return h_.Mode(); }

private:
	FactorMatrix<T> h_;
};

template<typename T>
using SPOM = LeafMatrix<T>;

typedef SPOM<complex<double>> SPOMcd;

typedef SPOM<double> SPOMd;

#endif //LEAFMATRIX_H
