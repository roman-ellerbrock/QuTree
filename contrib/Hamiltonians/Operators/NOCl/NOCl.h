#pragma once
#include "Core/stdafx.h"
#include "Core/Vector.h"
#include "TreeOperators/LeafFunction.h"
#include "TreeOperators/SumOfProductsOperator.h"

Tensorcd qv1(const LeafInterface& grid, const Tensorcd& Acoeffs);
Tensorcd qv2(const LeafInterface& grid, const Tensorcd& Acoeffs);
Tensorcd qv3(const LeafInterface& grid, const Tensorcd& Acoeffs);

Tensorcd qd1(const LeafInterface& grid, const Tensorcd& Acoeffs);
Tensorcd qd2(const LeafInterface& grid, const Tensorcd& Acoeffs);
Tensorcd qd3(const LeafInterface& grid, const Tensorcd& Acoeffs);
Tensorcd qd4(const LeafInterface& grid, const Tensorcd& Acoeffs);

Tensorcd qdminus(const LeafInterface& grid, const Tensorcd& Acoeffs);

Tensorcd ApplyNOpotential(const LeafInterface& grid, const Tensorcd& Acoeff);

Tensorcd ApplyW00(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW01(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW02(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW03(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW04(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW10(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW11(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW12(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW13(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW14(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW20(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW21(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW22(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW23(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW24(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW30(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW31(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW32(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW33(const LeafInterface& grid, const Tensorcd& Acoeff);
Tensorcd ApplyW34(const LeafInterface& grid, const Tensorcd& Acoeff);

class NOCl
	:public SOPcd
{
public:
	NOCl(const Tree& basis, bool potential);
	~NOCl() = default;

	Vectord getW(vector<double> coeffs, const LeafInterface & grid);

	void SpecialInitialize(const Tree & basis);
	void InitKin();
	void InitV(const Tree & basis);

protected:
	vector<Vectord> w;
	Vectord potno;
	bool potential;
};


