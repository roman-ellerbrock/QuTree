#include "SpinGroup.h"
#include <chrono>

void SpinGroup::Initialize(double par0, double par1, double par2, double par3) {
	last = false;
	alpha = par2 * 3.1415926538;
	if (par3 < 0) { last = true; }
}

Tensorcd SpinGroup::applyX(const Tensorcd& Acoeffs) const {
	return Acoeffs;
}

Tensorcd SpinGroup::ApplyX2(const Tensorcd& Acoeffs) const {
	return Acoeffs;
}

Tensorcd SpinGroup::ApplyP(const Tensorcd& Acoeffs) const {
	return Acoeffs;
}

Tensorcd SpinGroup::ApplyKin(const Tensorcd& Acoeffs) const {
	return Acoeffs;
}

void SpinGroup::InitSPF(Tensorcd& Acoeffs) const {
	const TensorDim& tdim = Acoeffs.Dim();
	size_t ntensor = tdim.getntensor();
	size_t dimpart = tdim.getdimpart();

	using namespace std::chrono;
	auto seed = system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	uniform_real_distribution<double> dist;

	if (dim == 2) {
		assert(ntensor == 2);
		assert(dimpart == 2);
		Acoeffs(0, 0) = cos(alpha);
		Acoeffs(1, 0) = sin(alpha);
		Acoeffs(0, 1) = cos(alpha);
		Acoeffs(1, 1) = -sin(alpha);
	} else {
		for (size_t n = 0; n < ntensor; ++n) {
			for (size_t i = 1; i < dimpart; ++i) {
				Acoeffs(i, n) = dist(gen);
			}
		}
		for (size_t i = 0; i < dimpart; ++i) {
			Acoeffs(i, 0) = 0.;
		}
		Acoeffs(0, 0) = 1.;

		if (last) {
			for (size_t i = 0; i < dimpart; ++i) {
				Acoeffs(i, 0) = 0.;
			}
			Acoeffs(dimpart - 1, 0) = 1.;
		}
	}

	GramSchmidt(Acoeffs);
}


