#include "CH3Potential.h"

extern "C"
{
// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
void __stdcall POTENTIALCH3(double* v, double* q, int* part);
#else
void potentialch3_(double *v, double *q, int *part);
#endif
}

double CH3Potential::evaluate(const Vectord& Xv, size_t part) const {
	double v = 0;
	Vectord Xv2(Xv);
	assert(Xv2.dim() == 12);
	int fpart = (int) part;

	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	POTENTIALCH3(&v, &Xv2(0), &fpart);
#else
	potentialch3_(&v, &Xv2(0), &fpart); //@TODO: did not compile
#endif
	// Cut-off for stability
	v = min(v, 5. / 27.2114);
	return v;
}
