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
		const TensorDim& dim = A.Dim();
		Tensorcd B(dim);
		InitSPF(B);
		for (size_t i = 0; i < dim.GetDimTot(); ++i) {
			A(i) = abs(B(i));
		}
	}

	virtual Tensorcd applyX(const Tensorcd& A)const = 0;
	virtual Tensorcd ApplyX2(const Tensorcd& A)const = 0;
	virtual Tensorcd ApplyP(const Tensorcd& A)const = 0;
	virtual Tensorcd ApplyKin(const Tensorcd& A)const = 0;
	virtual const Vectord& GetX()const = 0;
	virtual Vectord& GetX() = 0;
	virtual Tensorcd ToGrid(const Tensorcd& A)const = 0;
	virtual Tensorcd FromGrid(const Tensorcd& A)const = 0;
	virtual int oSQR()const = 0;
	virtual bool HasDVR()const = 0; // Tells wether a primitive basis does have a grid representation
};

