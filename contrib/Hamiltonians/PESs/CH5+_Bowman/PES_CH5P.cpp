#include "PES_CH5P.h"

extern "C"
{
	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	void __stdcall POTENTIALCH5P(double* v, double* q, int* part);
#endif
#ifdef __linux__
	void potentialch5p_(double* v, double* q, int* part);
#endif
}

PES_CH5P::PES_CH5P()
{
}

PES_CH5P::~PES_CH5P()
{
}

double PES_CH5P::Evaluate(const Vectord& Xv, size_t part)
{
	double v = 0;
	Vectord Xv2(Xv);
	assert(Xv2.Dim() == 12);
	int fortranpart = part;

	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	POTENTIALCH5P(&v, &Xv2(0), &fortranpart);
#endif
#ifdef __linux__
	potentialch5p_(&v, &Xv2(0), &fortranpart);
#endif
	// Cut-off for stability
	v = min(v, 5. / 27.2114);
	return v;
}
