#include "OCSHePotential.h"

extern "C"
{
#ifdef _MSC_VER
	void __stdcall INIT_POT();
#endif
#ifdef __linux__
		void init_pot_();
#endif
}

extern "C"
{
#ifdef _MSC_VER
  void __stdcall GET_POT_EN(double* x, double* y, double* z, double* ener);
#endif
#ifdef __linux__
  void get_pot_en_(double* x, double* y, double* z, double* ener);
#endif
}

OCSHePotential::OCSHePotential()
{
}

OCSHePotential::~OCSHePotential()
{
}

void OCSHePotential::Initialize()
{
#ifdef _MSC_VER
	INIT_POT();
#endif
#ifdef __linux__
	init_pot_();
#endif
}

double OCSHePotential::Evaluate(const Vectord& Xv, size_t part)
{
	double x = Xv(0);
	double y = Xv(1);
	double z = Xv(2);

	double ener;

#ifdef _MSC_VER
	GET_POT_EN(&x, &y, &z, &ener);
#endif
#ifdef __linux__
	get_pot_en_(&x, &y, &z, &ener);
#endif

	if(ener > 500.) ener = 500.;
	return (ener+70.)/cm;
}
