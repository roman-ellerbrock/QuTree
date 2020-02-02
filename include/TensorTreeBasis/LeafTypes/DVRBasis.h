#pragma once
#include "Core/Tensor.h"
#include "PrimitiveBasis.h"
#include "Core/Matrix.h"
#include "Core/Vector.h"
#include "Core/FactorMatrix.h"

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

	Vectord& GetX() override { return x_; }

	const Vectord& GetX() const override { return x_; }

	bool HasDVR() const override { return true; }

protected:

	// location operators
	Vectord x_;
	FactorMatrixd trafo_;

	// derivative operators
	FactorMatrixd kin_;
	FactorMatrixcd p_;

	// parameters that define the basis
	double omega_, r0_, wfr0_, wfomega_;
	int dim_;
};
