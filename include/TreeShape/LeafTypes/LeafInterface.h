#pragma once
#include"Core/Tensor.h"

class LeafInterface {
public:
	LeafInterface() = default;
	virtual ~LeafInterface() = default;

	virtual void initialize(double par0, double par1, double par2, double par3) = 0;
	virtual void initSPF(Tensorcd& A) const = 0;

	virtual void initSPF(Tensord& A) const {
		const TensorShape& dim = A.shape();
		Tensorcd B(dim);
		initSPF(B);
		for (size_t i = 0; i < dim.totalDimension(); ++i) {
			A(i) = abs(B(i));
		}
	}

	virtual void applyX(Tensorcd& xA, const Tensorcd& A) const = 0;
	virtual void applyX2(Tensorcd& x2A, const Tensorcd& A) const = 0;
	virtual void applyP(Tensorcd& pA, const Tensorcd& A) const = 0;
	virtual void applyKin(Tensorcd& kinA, const Tensorcd& A) const = 0;

	void identity(Tensorcd& IPhi, const Tensorcd& Phi) const {
		for (size_t i = 0; i < Phi.shape().totalDimension(); ++i) {
			IPhi(i) = Phi(i);
		}
	}

	virtual const Vectord& getX() const = 0;
	virtual Vectord& getX() = 0;
	virtual void toGrid(Tensorcd& UA, const Tensorcd& A) const = 0;
	virtual void fromGrid(Tensorcd& UA, const Tensorcd& A) const = 0;
	virtual int oSQR() const = 0;
	virtual bool hasDVR() const = 0; // Tells wether a primitive basis does have a grid representation
};

