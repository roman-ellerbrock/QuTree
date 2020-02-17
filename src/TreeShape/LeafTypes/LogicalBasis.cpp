#include "LogicalBasis.h"

// Logical Nodes do nothing
void LogicalBasis::Initialize(double par1, double par2, double par3, double par4) {
}

void LogicalBasis::InitSPF(Tensorcd& SPF)const {
	const TensorDim& tdim = SPF.Dim();
	size_t ntensor = tdim.GetNumTensor();
	size_t dimpart = tdim.LastBefore();
	SPF.Zero();
	assert(ntensor == dimpart);
	for (size_t n = 0; n < ntensor; n++) {
		SPF(n, n) = 1.;
	}
}

Tensorcd LogicalBasis::applyX(const Tensorcd& Acoeffs)const {
	assert(0);
	return Tensorcd();
}

Tensorcd LogicalBasis::ApplyX2(const Tensorcd& Acoeffs)const {
	assert(0);
	return Tensorcd();
}

Tensorcd LogicalBasis::ApplyP(const Tensorcd& Acoeffs)const {
	assert(0);
	return Tensorcd();
}

Tensorcd LogicalBasis::ApplyKin(const Tensorcd& Acoeffs)const {
	assert(0);
	return Tensorcd();
}


const Vectord& LogicalBasis::GetX()const {
	assert(0);
//	return Vectord();
}

Vectord& LogicalBasis::GetX(){
	assert(0);
//	return Vectord();
}

Tensorcd LogicalBasis::ToGrid(const Tensorcd& Acoeffs)const {
	assert(0);
	return Tensorcd();
}

Tensorcd LogicalBasis::FromGrid(const Tensorcd& Acoeffs)const {
	assert(0);
	return Tensorcd();
}

