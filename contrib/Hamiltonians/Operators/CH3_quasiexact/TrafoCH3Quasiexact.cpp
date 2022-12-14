#include "TrafoCH3Quasiexact.h"

extern "C" {
// @TODO: use cmake for this and ensure that it works for all relevant platforms
void inttocartch3qe_(double *q, double *x);
}

TrafoCH3Quasiexact::TrafoCH3Quasiexact(Vectord mass_)
		: mass(mass_) {}

Vectord TrafoCH3Quasiexact::transform(const Vectord& q) const {
	assert(q.dim() == 6);
	Vectord qq(q);
	Vectord Xv(12);
	inttocartch3qe_((double *) (&qq(0)), (double *) &Xv(0));
	return Xv;
}
