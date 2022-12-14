#include "PES_CH4SM.h"

extern "C"
{
	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	void __stdcall POTENTIALCH4SM(double* v, double* q, int* part);
	void __stdcall POTINICH4SM();
#endif
#ifdef __linux__
	void potentialch4sm_(double* v, double* q, int* part);
	void potinich4sm_();
#endif
}


void PES_CH4SM::Initialize()
{
#ifdef _MSC_VER
	POTINICH4SM();
#endif
#ifdef __linux__
	potinich4sm_();
#endif
}

PES_CH4SM::PES_CH4SM()
{
	Initialize();
}


PES_CH4SM::~PES_CH4SM()
{
}

double PES_CH4SM::Evaluate(const Vectord& Xv, size_t part)
{
	double v = 0;
	Vectord Xv2(Xv);
	assert(Xv2.Dim() == 9);
	int fortranpart = part;
	
	// @TODO: Use cmake and ensure it works on all relevant platforms
#ifdef _MSC_VER
	POTENTIALCH4SM(&v, &Xv2(0), &fortranpart);
#endif
#ifdef __linux__
	potentialch4sm_(&v, &Xv2(0), &fortranpart);
#endif
	// Cut-off for stability
	v = min(v, 5. / 27.2114);
	return v;
}
