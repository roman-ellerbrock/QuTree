#pragma once
#include "Tensor.h"
#include "PrimitiveBasis.h"
#include "Matrix.h"
#include "Vector.h"
#include "FactorMatrix.h"

class DVRBasis
	: public PrimitiveBasis {
public:
	DVRBasis() = default;
	explicit DVRBasis(int dim);
	virtual ~DVRBasis() override = default;

	void Initialize(double par0, double par1, double par2, double par3) override = 0;

	void InitSPF(Tensorcd& A) const override = 0;
	Tensorcd applyX(const Tensorcd& A) const override;
	Tensorcd ApplyX2(const Tensorcd& A) const override;
	Tensorcd ApplyP(const Tensorcd& A) const override;
	Tensorcd ApplyKin(const Tensorcd& A) const override;

	Tensorcd ToGrid(const Tensorcd& A) const override;
	Tensorcd FromGrid(const Tensorcd& A) const override;

	int oSQR() const override { return -1; }

	Vectord& GetX() override { return x; }

	const Vectord& GetX() const override { return x; }

	bool HasDVR() const override { return true; }

protected:

	// location operators
	Vectord x;
	FactorMatrixd trafo;

	// derivative operators
	FactorMatrixd kin;
	FactorMatrixcd p;

	// parameters that define the basis
	double omega, r0, wfr0, wfomega;
	int dim;
};
