#include "PES_CH4.h"


extern "C"
{
	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	void __stdcall POTENTIALCH4(double* v, double* q, int* part);
	void __stdcall POTINICH4();
#endif
#ifdef __linux__
	void potentialch4_(double* v, double* q, int* part);
	void potinich4_();
#endif
}


void PES_CH4::Initialize()
{
#ifdef _MSC_VER
	POTINICH4();
#endif
#ifdef __linux__
	potinich4_();
#endif
}

PES_CH4::PES_CH4()
{
	Initialize();
}


PES_CH4::~PES_CH4()
{
}

double PES_CH4::Evaluate(const Vectord& Xv, size_t part)
{
	double v = 0;
	Vectord Xv2(Xv);
	assert(Xv2.Dim() == 15);
	int fortranpart = part;

	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	POTENTIALCH4(&v, &Xv2(0), &fortranpart);
#endif
#ifdef __linux__
	potentialch4_(&v, &Xv2(0), &fortranpart);
#endif
	// Cut-off for stability
	v = min(v, 5. / 27.2114);
	return v;
}
