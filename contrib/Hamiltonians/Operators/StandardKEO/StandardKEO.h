#pragma once
#include "Hamiltonian.h"

class StandardKEO :
	public SumOfProductsOperator
{
public:
	StandardKEO(const mctdhBasis& basis, size_t dim);
	~StandardKEO();

	void SpecialInitialize(const mctdhBasis & basis, size_t dim);
};

