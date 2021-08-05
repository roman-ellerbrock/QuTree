#include "TreeShape/LeafTypes/SpinGroup.h"
#include <chrono>

void SpinGroup::initialize(double par0, double par1, double par2, double par3) {
    last_ = false;
    alpha_ = par2 * 3.1415926538;
	if (par3 < 0) { last_ = true; }
}

void SpinGroup::initSPF(Tensorcd& A) const {

	const TensorShape& shape = A.shape();
	if (shape.order() < 1 || shape[0] <= 1) {
		cerr << "Cannot initialize spin basis: primitive basis too small.\n";
		exit(1);
	}
	A(0, 0) = cos(alpha_);
	if (shape[0] > 1) {
		A(1, 0) = sin(alpha_);
	}
	if (A.shape().lastDimension() > 1) {
		A(0, 1) = -sin(alpha_);
		A(1, 1) = cos(alpha_);
	}

	if (A.shape().lastDimension() > 2 ) {
		using namespace std::chrono;
		uniform_real_distribution<double> dist;
		auto seed = system_clock::now().time_since_epoch().count();
		mt19937 gen(seed);
		size_t ntensor = A.shape().lastDimension();
		size_t dimpart = A.shape().lastBefore();
		for (size_t n = 2; n < ntensor; ++n) {
			for (size_t i = 2; i < dimpart; ++i) {
				A(i, n) = dist(gen);
			}
		}
	}

/*	using namespace std::chrono;
	uniform_real_distribution<double> dist;
	assert(ntensor == 2);
	assert(dimpart == 2);
	auto seed = system_clock::now().time_since_epoch().count();
	mt19937 gen(seed);
	size_t ntensor = tdim.lastDimension();
	size_t dimpart = tdim.lastBefore();
	if (dim_ == 2) {
		assert(ntensor == 2);
		assert(dimpart == 2);
		A(0, 0) = cos(alpha_);
		A(1, 0) = sin(alpha_);
		A(0, 1) = -sin(alpha_);
		A(1, 1) = cos(alpha_);
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

		if (last_) {
			for (size_t i = 0; i < dimpart; ++i) {
				A(i, 0) = 0.;
			}
			A(dimpart - 1, 0) = 1.;
		}
	}*/

	gramSchmidt(A);
}


