#include "hehepotential.h"

HeHePotential::HeHePotential()
{
}


HeHePotential::~HeHePotential()
{
}

extern "C"
{
#ifdef _MSC_VER
	void __stdcall HE_HE_POTENTIAL(double* r, double* e);
#endif
#ifdef __linux__
	void hehepotential_(double* r, double* e);
#endif
}

double HeHePotential::Evaluate(double x, size_t part)
{
	double ener;
#ifdef _MSC_VER
	HEHEPOTENTIAL(&x, &ener);
#endif
#ifdef __linux__
	hehepotential_(&x, &ener);
#endif
	return ener;
}
