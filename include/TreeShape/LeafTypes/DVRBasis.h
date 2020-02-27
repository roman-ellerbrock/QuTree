#pragma once
#include "Core/Tensor.h"
#include "LeafInterface.h"
#include "Core/Matrix.h"
#include "Core/Vector.h"

class DVRBasis
	: public LeafInterface {
public:
	DVRBasis() = default;
	explicit DVRBasis(int dim);
	virtual ~DVRBasis() override = default;

	void Initialize(double par0, double par1, double par2, double par3) override = 0;

	void InitSPF(Tensorcd& A) const override = 0;
	void  applyX(Tensorcd& uA, const Tensorcd& A) const override;
	void  applyX2(Tensorcd& uA, const Tensorcd& A) const override;
	void  applyP(Tensorcd& uA, const Tensorcd& A) const override;
	void  applyKin(Tensorcd& uA, const Tensorcd& A) const override;

	void ToGrid(Tensorcd& uA, const Tensorcd& A) const override;
	void FromGrid(Tensorcd& uA, const Tensorcd& A) const override;

	int oSQR() const override { return -1; }

	Vectord& GetX() override { return x_; }

	const Vectord& GetX() const override { return x_; }

	bool HasDVR() const override { return true; }

protected:

	// location operators
	Vectord x_;
	Matrixd trafo_;

	// derivative operators
	Matrixd kin_;
	Matrixcd p_;

	// parameters that define the basis
	double omega_, r0_, wfr0_, wfomega_;
	int dim_;
};
