#include "CDVRModelV2.h"



CDVRModelV2::CDVRModelV2()
{
}


CDVRModelV2::~CDVRModelV2()
{
}

double CDVRModelV2::Evaluate(const Vectord & Xv_in, size_t part)
{
	constexpr double cm = 219474.6313705;
	constexpr double lambda = 50. / cm;
	constexpr double omega = 500. / cm;
	constexpr double omegamorse = 3000. / cm;
	// DOF
	size_t f = 8;
	double v = 0;
	double D = 6 * omegamorse;
	double alpha = sqrt(omegamorse * omegamorse / (2. * D));

	Vectord vmorse(2);
	Vectord Xv(Xv_in);
	Xv(0)   = sqrt(0.95) * Xv_in(0) + sqrt(0.05) * Xv_in(f-1);
	Xv(f-1) = sqrt(0.05) * Xv_in(0) - sqrt(0.95) * Xv_in(f-1);
	Xv(1)   = sqrt(0.95) * Xv_in(1) + sqrt(0.05) * Xv_in(f-2);
	Xv(f-2) = sqrt(0.05) * Xv_in(1) - sqrt(0.95) * Xv_in(f-2);

	// Seperable parts
	for (int i = 0; i < 2; ++i)
	{
		vmorse(i) = D * pow(1. - exp(-alpha*Xv(i)), 2);
//		vmorse(i) = 0.5*pow(omegamorse*Xv(i), 2);
		v += vmorse(i);
	}

	// "Bend"-coordinates 2-f coupled in a ring
	for (int i = 2; i < f; i++)
		v += 0.5*omega*omega*Xv(i)*Xv(i);

	// Couplings

//	for (int k = 0; k < 2; ++k)
//		for (int i = 2; i < f; ++i)
//			v += 0.3*omega*omega*Xv(i)*Xv(k); 

	for (int i = 2; i < f - 1; ++i)
		v += lambda*lambda*Xv(i)*Xv(i+1); 
	v += lambda*lambda*Xv(f - 1)*Xv(2);

	return v;
}

