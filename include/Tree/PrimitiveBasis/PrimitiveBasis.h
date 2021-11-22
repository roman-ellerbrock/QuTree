#pragma once
#include "Tensor/Tensor.h"

struct BasisParameters {
	/**
	 * \class Parameter container of a leaf.
	 */
	~BasisParameters() = default;

	void read(istream& file) {
		file >> par0_;
		file >> par1_;
		file >> par2_;
		file >> par3_;
	}

	[[nodiscard]] double omega() const { return par0_; }
	[[nodiscard]] double r0() const { return par1_; }
	[[nodiscard]] double wfomega() const { return par2_; }
	[[nodiscard]] double wfr0() const { return par3_; }

	double par0_, par1_, par2_, par3_;
};

class LeafInterface {
public:
	LeafInterface() = default;
	~LeafInterface() = default;

	virtual void resize(size_t dim) = 0;
	virtual void initialize(size_t dim, const BasisParameters& par) = 0;
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

	virtual const Tensord& getX() const = 0;
	virtual Tensord& getX() = 0;

	virtual void toGrid(Tensorcd& UA, const Tensorcd& A) const = 0;
	virtual void fromGrid(Tensorcd& UA, const Tensorcd& A) const = 0;

	virtual bool hasDVR() const = 0; // Tells wether a primitive basis does have a grid representation

	BasisParameters par_;
};

