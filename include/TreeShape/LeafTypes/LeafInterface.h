#pragma once
#include"Core/Tensor.h"

class LeafInterface
{
public:
	LeafInterface() = default;
	virtual ~LeafInterface() = default;

	virtual void Initialize(double par0, double par1, double par2, double par3) = 0;
	virtual void InitSPF(Tensorcd& A)const = 0;

	virtual void InitSPF(Tensord& A)const {
		const TensorShape& dim = A.shape();
		Tensorcd B(dim);
		InitSPF(B);
		for (size_t i = 0; i < dim.totalDimension(); ++i) {
			A(i) = abs(B(i));
		}
	}

	virtual void applyX(Tensorcd& xA, const Tensorcd& A)const = 0;
	virtual void applyX2(Tensorcd& x2A, const Tensorcd& A)const = 0;
	virtual void applyP(Tensorcd& pA, const Tensorcd& A)const = 0;
	virtual void applyKin(Tensorcd& kinA, const Tensorcd& A)const = 0;
	virtual const Vectord& GetX()const = 0;
	virtual Vectord& GetX() = 0;
	virtual void ToGrid(Tensorcd& UA, const Tensorcd& A)const = 0;
	virtual void FromGrid(Tensorcd& UA, const Tensorcd& A)const = 0;
	virtual int oSQR()const = 0;
	virtual bool HasDVR()const = 0; // Tells wether a primitive basis does have a grid representation
};

