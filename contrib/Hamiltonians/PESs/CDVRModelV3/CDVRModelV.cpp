#include "CDVRModelV.h"



CDVRModelV3::CDVRModelV3()
	:simple(false)
{
}


CDVRModelV3::~CDVRModelV3()
{
}

vector<MultiParticleOperator> CDVRModelV3::Build(
	const	mctdhBasis& basis)
{
	cout << "=== Potential ===" << endl;
	cout << "Initializing two adiabatic potentials." << endl;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::applyX;

	vector<MultiParticleOperator> ms;
	size_t f = basis.nLeaves();

	PotentialOperator v1(3, 0);
	MultiParticleOperator M1;
	M1.SetV(v1);
	ms.push_back(M1);

	if (!simple) {
		PotentialOperator v2(3, 1);
		MultiParticleOperator M2;
		M2.push_back(x, f - 1);
		M2.SetV(v2);
		ms.push_back(M2);
	}

	cout << "Added " << ms.size() << " parts to the Potential." << endl;
	cout << "=== Leaving Potential ===" << endl;
	return ms;
}

double CDVRModelV3::Evaluate(const Vectord & Xv, size_t part)
{
	constexpr double cm = 219474.6313705;
	constexpr double lambda = 2000. / cm;
	constexpr double omega = 4000. / cm;
	double v = 0;

	if (part == 0)
	{
		for (int i = 0; i < 3; i++)
			v += 0.5*omega*omega*Xv(i)*Xv(i);
	
		v += lambda*lambda*Xv(0)*Xv(1);
		v += lambda*lambda*Xv(1)*Xv(2);
	}
	else if (part == 1)
	{
		// x_2 * x_3
		v += lambda*lambda*(Xv(2) + Xv(0));
	}
	else
	{
		cerr << "Wrong part number." << endl;
		exit(1);
	}

	return v;
}

