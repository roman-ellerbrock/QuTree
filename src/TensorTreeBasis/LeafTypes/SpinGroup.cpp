#include "SpinGroup.h"
#include <chrono>

void SpinGroup::Initialize(double par0, double par1, double par2, double par3) {
	last = false;
	alpha = par2 * 3.1415926538;
	if (par3 < 0) { last = true; }
}

Tensorcd SpinGroup::applyX(const Tensorcd& A) const {
	return A;
}

Tensorcd SpinGroup::ApplyX2(const Tensorcd& A) const {
	return A;
}

Tensorcd SpinGroup::ApplyP(const Tensorcd& A) const {
	return A;
}

Tensorcd SpinGroup::ApplyKin(const Tensorcd& A) const {
	return A;
}

void SpinGroup::InitSPF(Tensorcd& A) const {
	const TensorDim& tdim = A.Dim();
	size_t ntensor = tdim.GetNumTensor();
	size_t dimpart = tdim.GetDimPart();

	using namespace std::chrono;
	auto seed = system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	uniform_real_distribution<double> dist;

	if (dim == 2) {
		assert(ntensor == 2);
		assert(dimpart == 2);
		A(0, 0) = cos(alpha);
		A(1, 0) = sin(alpha);
		A(0, 1) = -sin(alpha);
		A(1, 1) = cos(alpha);
	} else {
		for (size_t n = 0; n < ntensor; ++n) {
			for (size_t i = 1; i < dimpart; ++i) {
				A(i, n) = dist(gen);
			}
		}
		for (size_t i = 0; i < dimpart; ++i) {
			A(i, 0) = 0.;
		}
		A(0, 0) = 1.;

		if (last) {
			for (size_t i = 0; i < dimpart; ++i) {
				A(i, 0) = 0.;
			}
			A(dimpart - 1, 0) = 1.;
		}
	}

	GramSchmidt(A);
}


