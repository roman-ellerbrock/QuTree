#include "Tree/PrimitiveBasis/SpinGroup.h"
#include "Util/QMConstants.h"
#include "Tensor/Tensor"

namespace JordanWigner {
	Tensorcd identity() {
		return identitycd({2, 2});
	}

	Tensorcd sigmaX() {
		Tensorcd x({2, 2});
		x(1, 0) = 1.;
		x(0, 1) = 1.;
		return x;
	}

	Tensorcd sigmaY() {
		Tensorcd y({2, 2});
		y(1, 0) = QM::im;
		y(0, 1) = -QM::im;
		return y;
	}

	Tensorcd sigmaZ() {
		Tensorcd z({2, 2});
		z(0, 0) = 1.;
		z(1, 1) = -1.;
		return z;
	}

	Tensorcd sigmaPlus() {
		/// Ref. [1] Eq. (19)
		Tensorcd s({2, 2});
		s(1, 0) = 1.;
		return s;
//		return 0.5 * (sigmaX() - QM::im * sigmaY());
	}

	Tensorcd sigmaMinus() {
		/// Ref. [1] Eq. (19)
		Tensorcd s({2, 2});
		s(0, 1) = 1.;
		return s;
//		return 0.5 * (sigmaX() + QM::im * sigmaY());
	}
}

void SpinGroup::initialize(const BasisParameters& par) {
	par_ = par;
	alpha_ = par_.par2_ * QM::pi;
}

void SpinGroup::occupy(Tensorcd& A) const {

	const TensorShape& shape = A.shape_;
	if (shape.order() < 1 || shape[0] <= 1) {
		cerr << "Cannot initialize spin basis: primitive basis too small.\n";
		exit(1);
	}

	A(0, 0) = cos(alpha_);

	if (shape[0] > 1) {
		A(1, 0) = sin(alpha_);
	}

	if (A.shape_.lastDimension() > 1) {
		A(0, 1) = -sin(alpha_);
		A(1, 1) = cos(alpha_);
	}

	if (A.shape_.lastDimension() > 2 ) {
		using namespace std::chrono;
		uniform_real_distribution<double> dist;
		auto seed = system_clock::now().time_since_epoch().count();
		mt19937 gen(seed);
		size_t ntensor = A.shape_.lastDimension();
		size_t dimpart = A.shape_.lastBefore();
		for (size_t n = 2; n < ntensor; ++n) {
			for (size_t i = 2; i < dimpart; ++i) {
				A(i, n) = dist(gen);
			}
		}
	}

	gramSchmidt(A);
}


