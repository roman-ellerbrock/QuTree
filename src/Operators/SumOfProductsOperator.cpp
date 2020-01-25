#include "SumOfProductsOperator.h"

SumOfProductsOperator multAB(const SumOfProductsOperator& A,
	const SumOfProductsOperator& B)
{
	SumOfProductsOperator C;

	for (size_t i = 0; i < A.size(); i++)
	{
		const MultiParticleOperator Ai = A(i);
		auto acoeff = A.Coeff(i);
		for (size_t j = 0; j < B.size(); j++)
		{
			const MultiParticleOperator Bj = B(j);
			auto bcoeff = B.Coeff(j);
			auto ccoeff = acoeff*bcoeff;
			C.push_back(Ai*Bj, ccoeff);
		}
	}

	return C;
}

