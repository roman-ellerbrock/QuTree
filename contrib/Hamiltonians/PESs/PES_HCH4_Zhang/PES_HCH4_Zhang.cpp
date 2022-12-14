#include "PES_HCH4_Zhang.h"

extern "C"
{
	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	void __stdcall POTENTIALHCH4ZHANG(double* v, double* q, int* part);
	void __stdcall POTINIHCH4ZHANG();
#endif
#ifdef __linux__
	void potentialhch4zhang_(double* v, double* q, int* part);
	void potinihch4zhang_();
#endif
}

PES_HCH4_Zhang::PES_HCH4_Zhang()
{
	Initialize();
}

PES_HCH4_Zhang::~PES_HCH4_Zhang()
{
}

void PES_HCH4_Zhang::Initialize()
{
#ifdef _MSC_VER
	POTINIHCH4ZHANG();
#endif
#ifdef __linux__
	potinihch4zhang_();
#endif
}

double PES_HCH4_Zhang::Evaluate(const Vectord& Xv, size_t part)
{
	double v = 0;
	Vectord Xv2(Xv);
	assert(Xv2.Dim() == 18);
	int fortranpart = part;

	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	POTENTIALHCH4ZHANG(&v, &Xv2(0), &fortranpart);
#endif
#ifdef __linux__
	potentialhch4zhang_(&v, &Xv2(0), &fortranpart);
#endif
	// Cut-off for stability
	v = min(v, 5. / 27.2114);
	return v;
}
