//
// Created by Roman Ellerbrock on 8/9/20.
//

#ifndef SPINBOSON_H
#define SPINBOSON_H
#include "TreeOperators/LeafFunction.h"
#include "TreeOperators/SumOfProductsOperator.h"

namespace Operator {

	SOPcd Exciton(const string& filename, const Tree& tree);

	class Excite: public LeafOperatorcd {
	public:
		Excite(size_t from, size_t to)
			: from_(from), to_(to) {}
		~Excite() = default;

		void apply(const LeafInterface& grid, Tensorcd& PPhi,
			const Tensorcd& Phi) const override {

			const TensorShape& tdim = Phi.shape();
			PPhi.zero();
			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				PPhi(to_, n) = Phi(from_, n);
			}
		}

	private:
		size_t from_;
		size_t to_;
	};

	class VectorExcite: public LeafOperatorcd {
	public:
		VectorExcite(size_t from, const Vectord& a)
			: a_(a), from_(from) {}
		~VectorExcite() {}

		void apply(const LeafInterface& grid, Tensorcd& PPhi,
			const Tensorcd& Phi) const override {

			const TensorShape& tdim = Phi.shape();
			PPhi.zero();
			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				for (size_t i = 0; i < tdim.lastBefore(); ++i) {
					PPhi(i, n) = a_(i) * Phi(from_, n);
				}
			}
		}

	private:
		Vectord a_;
		size_t from_;
	};

	class ExtendedX: public LeafOperatorcd {
	public:
		ExtendedX(size_t idx1, size_t idx2)
			: idx1_(idx1), idx2_(idx2) {}
		~ExtendedX() = default;

		void apply(const LeafInterface& grid, Tensorcd& PPhi,
			const Tensorcd& Phi) const override {

			const TensorShape& tdim = Phi.shape();
			PPhi.zero();
			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				PPhi(idx1_, n) = Phi(idx2_, n);
				PPhi(idx2_, n) = Phi(idx1_, n);
			}
		}

	private:
		size_t idx1_;
		size_t idx2_;
	};

	class Projector: public LeafOperatorcd {
	public:
		Projector(size_t idx)
			: idx_(idx) {}
		~Projector() = default;

		void apply(const LeafInterface& grid, Tensorcd& PPhi,
			const Tensorcd& Phi) const override {

			const TensorShape& tdim = Phi.shape();
			PPhi.zero();
			assert(idx_ < tdim.lastBefore());
			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				PPhi(idx_, n) = Phi(idx_, n);
			}
		}

	private:
		size_t idx_;
	};

	class VectorProjector: public LeafOperatorcd {
	public:
		explicit VectorProjector(const Vectord& U) : U_(U) {
		}

		void apply(const LeafInterface& grid, Tensorcd& PPhi,
			const Tensorcd& Phi) const override {
			const TensorShape& tdim = Phi.shape();

			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				for (size_t i = 0; i < tdim.lastBefore(); ++i) {
					PPhi(i, n) = U_(i) * Phi(i, n);
				}
			}
		}

	protected:
		Vectord U_;
	};

}

#endif //SPINBOSON_H
