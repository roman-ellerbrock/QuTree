#pragma once
#include "Tensor/Tensor.h"
#include "Tree/PrimitiveBasis/BasisParameters.h"

class PrimitiveBasis {
public:
	PrimitiveBasis() : par_() {}
	virtual ~PrimitiveBasis() = default;

	virtual void initialize(const BasisParameters& par) = 0;
	virtual void occupy(Tensorcd& A) const = 0;

	virtual void occupy(Tensord& A) const {
		const TensorShape& dim = A.shape_;
		Tensorcd B(dim);
		occupy(B);
		for (size_t i = 0; i < dim.totalDimension(); ++i) {
			A(i) = abs(B(i));
		}
	}

	virtual void applyX(Tensorcd& xA, const Tensorcd& A) const = 0;
	virtual void applyX2(Tensorcd& x2A, const Tensorcd& A) const = 0;
	virtual void applyP(Tensorcd& pA, const Tensorcd& A) const = 0;
	virtual void applyKin(Tensorcd& kinA, const Tensorcd& A) const = 0;

	void identity(Tensorcd& IPhi, const Tensorcd& Phi) const {
		for (size_t i = 0; i < Phi.shape_.totalDimension(); ++i) {
			IPhi(i) = Phi(i);
		}
	}

	void write(ostream& os) const {
		os << par_.dim_ << "\t" << par_.type_ << "\t" << par_.mode_ << "\n";
	}

	virtual const Tensord& getX() const = 0;
	virtual Tensord& getX() = 0;

	virtual void toGrid(Tensorcd& UA, const Tensorcd& A) const = 0;
	virtual void fromGrid(Tensorcd& UA, const Tensorcd& A) const = 0;

	virtual bool hasDVR() const = 0; // Tells wether a primitive basis does have a grid representation

	BasisParameters par_;
};

