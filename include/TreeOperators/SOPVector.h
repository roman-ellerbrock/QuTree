#pragma once
#include "SumOfProductsOperator.h"


template<typename T>
class SOPVector: public std::vector<SOP<T>>
/**
 * \class SOPVector
 * \ingroup Operators
 * \brief Dressed up vector of SOPs
 */
{
public:
	using vector<SOP<T>>::push_back;
	using vector<SOP<T>>::insert;
	void append(const SOP<T>& A) {
		this->push_back(A);
	}

	void append(const SOPVector<T>& A) {
		insert(this->end(), A.begin(), A.end());
	}

	//////////////////////////////////////////////////////////////////////
	/// Operators
	//////////////////////////////////////////////////////////////////////
	/// multiply with coefficient
	friend SOPVector<T> operator*(const MLO<T>& M,
		const SOPVector<T>& A) {
		SOPVector C;
		for (size_t i = 0; i < A.size(); i++) {
			C.push_back(M * A[i]);
		}
		return C;
	}
};

typedef SOPVector<complex<double>> SOPVectorcd;