#include "CDVRModelV.h"

CDVRModelV::CDVRModelV(size_t f_, bool coupling)
:f(f_), coupling_(coupling)
{}

double CDVRModelV::evaluate(const Vectord & Xv, size_t part)const
{
	constexpr double cm = 219474.6313705;
	constexpr double lambda = 2000. / cm;
	constexpr double omega = 4000. / cm;
	double v = 0;
	for (int i = 0; i < f; i++)
		v += 0.5*omega*omega*Xv(i)*Xv(i);

/*	double alpha = 500. / cm;
	for (int i = 0; i < f; i++)
		v += 0.5*alpha*alpha*Xv(i)*Xv(i)*Xv(i)*Xv(i);
*/
	if (coupling_) {
		for (size_t i = 0; i < f; ++i) {
			v += lambda * lambda * Xv(i) * Xv((i + 1) % f);
		}
	}

	return v;
}
