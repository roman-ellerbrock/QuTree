#include "Model1.h"



Model1::Model1()
	:simple(false)
{
}


Model1::~Model1()
{
}

vector<MultiParticleOperator> Model1::Build(
	const	mctdhBasis& basis)
{
	cout << "=== Potential ===" << endl;
	cout << "Initializing one adiabatic potential(s)." << endl;
	function<Tensorcd(const PrimitiveBasis&, const Tensorcd&)> x = &PrimitiveBasis::applyX;

	vector<MultiParticleOperator> ms;
	size_t f = basis.nLeaves();

	PotentialOperator v1(f, 0);
	MultiParticleOperator M1;
	M1.SetV(v1);
	ms.push_back(M1);

	cout << "Added " << ms.size() << " parts to the Potential." << endl;
	cout << "=== Leaving Potential ===" << endl;
	return ms;
}

double distort2d(double x, double y)
{
	constexpr double cm = 219474.6313705;
	double pi = 3.14159265359;
	double spatial_f = 0.06*2*pi;
	double amp = 100./cm;
	double distortion = amp * sin(spatial_f*x)*sin(spatial_f*y);
	return distortion;

}

double distort(double x)
{
	constexpr double cm = 219474.6313705;
	double pi = 3.14159265359;
	double spatial_f = 0.065*2*pi;
	double amp = 120./cm;
	double distortion = amp * sin(pow(spatial_f*x,1));
	return distortion;

}

double Model1::Evaluate(const Vectord & Xv, size_t part)
{
	constexpr double cm = 219474.6313705;
	constexpr size_t f = 4;
	double v = 0;
	constexpr size_t switcher = 1;

	Vectord omegas(4);
	omegas(0)=2800./cm;
	omegas(1)=1500./cm;
	omegas(2)=2800./cm;
	omegas(3)=1500./cm;

	if (switcher == 0) {
	
		Vectord y(4);
		double mix = 0.05;
		y(0)=sqrt(1-mix)*Xv(0)+sqrt(mix)*Xv(1);
		y(1)=-sqrt(mix)*Xv(0)+sqrt(1-mix)*Xv(0);
	
		y(2)=sqrt(1-mix)*Xv(2)+sqrt(mix)*Xv(3);
		y(3)=-sqrt(mix)*Xv(2)+sqrt(1.-mix)*Xv(3);

		constexpr double height = 10000 / cm;
		for (int k = 0; k < f; k++)
		{
			double freq = omegas(k) * omegas(k);
			double alpha = 0.5 * freq / height;
			v +=  height * (1. - exp( -alpha*y(k)*y(k)));
		}
	}else if (switcher == 1){
		Vectord y(4);
		double mix = 0.00;
		y(0)=sqrt(1-mix)*Xv(0)+sqrt(mix)*Xv(1);
		y(1)=Xv(1);
	
		y(2)=sqrt(1-mix)*Xv(2)+sqrt(mix)*Xv(3);
		y(3)=Xv(3);
		for (size_t k = 0; k < f; k++) {
			v += 0.5 * omegas(k)*omegas(k)*y(k)*y(k);
		}
//		v += distort2d(y(0), y(1));
//		v += distort2d(y(2), y(3));
//
		v += distort((y(0)+y(1))/sqrt(2.));
		v += distort((y(2)+y(3))/sqrt(2.));
		
//		v += distort((y(0)+y(1)+y(2)+y(3))/2.);
	}else{
		cerr << "PES: fatal error." << endl;
		exit(1);
	}

	return v;
}

