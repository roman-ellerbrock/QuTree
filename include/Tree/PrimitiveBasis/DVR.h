#pragma once
#include "Tensor/Tensor.h"
#include "PrimitiveBasis.h"

class DVR
	: public LeafInterface {
public:
	DVR() = default;
	~DVR() = default;

	void resize(size_t dim) override;
	void initialize(size_t dim, const BasisParameters& par) override;

	void occupy(Tensorcd& A) const override;

	void applyX(Tensorcd& uA, const Tensorcd& A) const override;
	void applyX2(Tensorcd& uA, const Tensorcd& A) const override;
	void applyP(Tensorcd& uA, const Tensorcd& A) const override;
	void applyKin(Tensorcd& uA, const Tensorcd& A) const override;

	void toGrid(Tensorcd& uA, const Tensorcd& A) const override;
	void fromGrid(Tensorcd& uA, const Tensorcd& A) const override;

	Tensord& getX() override { return x_; }

	[[nodiscard]] const Tensord& getX() const override { return x_; }

	[[nodiscard]] bool hasDVR() const override { return true; }


	[[nodiscard]] virtual double transformX(double x, bool forth) const {
		return x;
	}

	void shift(Tensord& x, double delta) const {
		for (size_t i = 0; i < x.shape_.totalDimension(); ++i) {
			x(i) = transformX(x(i), true);
			x(i) += delta;
			x(i) = transformX(x(i), false);
		}
	}

	Tensord x_;
	Tensorcd trafo_;
	Tensorcd w_;

	Tensorcd kin_;
	Tensorcd p_;

protected:
	[[nodiscard]] virtual Tensorcd buildX(size_t dim) const {
		return Tensorcd({dim, dim});
	}

	[[nodiscard]] virtual Tensorcd buildP(size_t dim) const {
		return Tensorcd({dim, dim});
	}

	[[nodiscard]] virtual Tensorcd buildKin(size_t dim) const {
		return Tensorcd({dim, dim});
	}

	[[nodiscard]] virtual Tensorcd buildW(size_t dim) const {
		return Tensorcd({dim});
	}
};
