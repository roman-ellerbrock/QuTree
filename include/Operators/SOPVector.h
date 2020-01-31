#pragma once
#include "SumOfProductsOperator.h"


class SOPVector: public std::vector<SOP>
/**
 * \class This class represents a vector of SOPs.
 * \ingroup Operators
 * \brief Dressed up vector of SOPs
 */
{
public:
	void append(const SOP& A) {
		push_back(A);
	}

	void append(const SOPVector& A) {
		insert(end(), A.begin(), A.end());
	}

	//////////////////////////////////////////////////////////////////////
	// Operators
	//////////////////////////////////////////////////////////////////////
	// multiply with coefficient
	friend SOPVector operator*(const MPO& M,
		const SOPVector& A) {
		SOPVector C;
		for (size_t i = 0; i < A.size(); i++) {
			C.push_back(M * A[i]);
		}
		return C;
	}
};

typedef SOPVector sopList;