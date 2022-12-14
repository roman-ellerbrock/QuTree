#include "NOClPotential.h"

NOClPotential::NOClPotential(bool dissociation_)
	: dissociation(dissociation_) {}

double NOClPotential::evaluate(const Vectord& Xv, size_t part) const {
	if (dissociation) {
		return evaluateS1(Xv, part);
	} else {
		return evaluateGS(Xv, part);
	}
}

double NOClPotential::evaluateGS(const Vectord& Xv, size_t part) const {
	const double mh = 1836.;
//	const double mN = 14.0 * mh;
	const double mCl = 35.403 * mh;
	const double mO = 16.0 * mh;

//	const double mud = 1. / (mN + mO) + 1. / mCl;
//	const double muv = 1. / mN + 1. / mO;
	const double mN = 25529.;
	const double muv = 1./13615.5;
	const double mud = 1./29456.5;
//	cout << "muv: " << 1./muv << endl;
//	cout << "mud: " << 1./mud << endl;
//	cout << "muN: " << mN << endl;
//	getchar();

	const double smuv = sqrt(1. / muv);
	const double smud = sqrt(1. / mud);

	const double mfracO = mO / (mO + mN);
	const double mfracN = mO / (mO + mN);

	double rv = Xv(0) / smuv;
	double rd = Xv(1) / smud;
	double cost = Xv(2);

	Vectord r(3);
	r(0) = rv;
	r(1) = sqrt(pow(mfracO * rv, 2) + rd * rd + 2 * rv * rd * cost * mfracO);
	r(2) = sqrt(pow(mfracN * rv, 2) + rd * rd - 2 * rv * rd * cost * mfracN);

	// Shift r by equilibrium values
	Vectord re(3);
	re(0) = 2.155;
	re(1) = 3.729;
	re(2) = 4.989;
	r = r - re;

	Vectord k(3);
	k(0) = 0.8987;
	k(1) = 0.0874;
	k(2) = 0.1137;
	double k23 = 0.0122;

	double V = 0;
	for (size_t i = 0; i < 3; ++i) {
		V += 0.5 * k(i) * pow(r(i), 2);
	}
	V += k23 * r(1) * r(2);
	V = min(5. / 27.2114, V);
	return V;
}

double NOClPotential::evaluateS1(const Vectord& Xv, size_t part) const {

	// mu stuff
//	const double mh = 1836.;
//	const double mN = 14.0 * mh;
//	const double mCl = 35.403 * mh;
//	const double mO = 16.0 * mh;
//	const double mud = 1. / (mN + mO) + 1. / mCl;
//	const double smud = sqrt(1. / mud);
	const double mN = 25529.;
	const double muv = 1./13615.5;
	const double mud = 1./29456.5;
	const double smud = sqrt(1. / mud);

	// Dissociation coordinate
	double xd = Xv(1) / smud;
	double rd = 1. - exp(-1.5 * (xd - 4.315));
/*	if (rd < (1. - exp(-1.5*(3.5 - 4.315))))
	{
		rd = (1. - exp(-1.5*(3.5 - 4.315)));
	}
	*/

	// Vibration coordinate
//	const double muv = 1. / mN + 1. / mO;
	const double smuv = sqrt(1. / muv);

	// Dissociation coordinate
	const double rv = Xv(0) / smuv - 2.136;

	const double pi = 3.14159265359;
	const double a = exp(-1.1 * Xv(2)) - exp(-1.1 * cos(pi * 127.4 / 180));

	double v =
		0.6816 * pow(rv, 2) - 0.9123 * pow(rv, 3) + 0.4115 * pow(rv, 4) +
			((0.0384816 + 0.0247875 * a + 0.0270933 * pow(a, 2) + 0.00126791 * pow(a, 3) +
				0.00541285 * pow(a, 4) + 0.0313629 * pow(a, 5) + 0.0172449 * pow(a, 6)) +
				(0.00834237 + 0.00398713 * a + 0.00783319 * pow(a, 2) + 0.0294887 * pow(a, 3) -
					0.0154387 * pow(a, 4) - 0.0621984 * pow(a, 5) - 0.0337951 * pow(a, 6)) * rd +
				(0.00161625 - 0.000156328 * a - 0.0189982 * pow(a, 2) - 0.00753297 * pow(a, 3) +
					0.00383665 * pow(a, 4) - 0.00758225 * pow(a, 5) - 0.00904493 * pow(a, 6)) * pow(rd, 2) +
				(-0.00101009 + 0.000619148 * a - 0.0149812 * pow(a, 2) - 0.0199722 * pow(a, 3) +
					0.00873013 * pow(a, 4) + 0.0376118 * pow(a, 5) + 0.0221523 * pow(a, 6)) * pow(rd, 3) +
				(-0.000368904 + 0.000164037 * a - 0.00331809 * pow(a, 2) - 0.00567787 * pow(a, 3) +
					0.00268662 * pow(a, 4) + 0.0134483 * pow(a, 5) + 0.0084585 * pow(a, 6)) * pow(rd, 4) +
				(-0.0558666 - 0.0276576 * a + 0.0934932 * pow(a, 2) - 0.0295638 * pow(a, 3) -
					0.15436 * pow(a, 4) + 0.0796119 * pow(a, 5) + 0.135121 * pow(a, 6) +
					(0.0582169 + 0.0384404 * a + 0.078114 * pow(a, 2) + 0.185556 * pow(a, 3) -
						0.0641656 * pow(a, 4) - 0.175976 * pow(a, 5) - 0.0104994 * pow(a, 6)) * rd +
					(0.052285 + 0.0472724 * a - 0.216008 * pow(a, 2) - 0.147775 * pow(a, 3) +
						0.349283 * pow(a, 4) + 0.28458 * pow(a, 5) + 0.00384449 * pow(a, 6)) * pow(rd, 2) +
					(0.0212609 + 0.0290597 * a - 0.109124 * pow(a, 2) + 0.0310445 * pow(a, 3) +
						0.262513 * pow(a, 4) - 0.250653 * pow(a, 5) - 0.369466 * pow(a, 6)) * pow(rd, 3) +
					(0.00334178 + 0.0039061 * a - 0.0110452 * pow(a, 2) + 0.0582029 * pow(a, 3) +
						0.0679524 * pow(a, 4) - 0.16459 * pow(a, 5) - 0.165337 * pow(a, 6)) * pow(rd, 4)) * rv +
				(-0.163186 - 0.180535 * a + 0.04692 * pow(a, 2) + 0.471673 * pow(a, 3) +
					0.403267 * pow(a, 4) - 0.718071 * pow(a, 5) - 0.761199 * pow(a, 6) +
					(-0.0290674 - 0.0136172 * a - 0.108952 * pow(a, 2) - 1.68269 * pow(a, 3) -
						1.2673 * pow(a, 4) + 3.17648 * pow(a, 5) + 2.92793 * pow(a, 6)) * rd +
					(0.121228 + 0.202308 * a + 0.483613 * pow(a, 2) + 1.29095 * pow(a, 3) -
						0.174483 * pow(a, 4) - 2.4605 * pow(a, 5) - 1.36597 * pow(a, 6)) * pow(rd, 2) +
					(0.107233 + 0.115213 * a - 0.366102 * pow(a, 2) + 0.812662 * pow(a, 3) +
						1.76038 * pow(a, 4) - 1.19665 * pow(a, 5) - 1.77392 * pow(a, 6)) * pow(rd, 3) +
					(0.0232767 + 0.0304932 * a - 0.19455 * pow(a, 2) - 0.0307517 * pow(a, 3) +
						0.539365 * pow(a, 4) + 0.120203 * pow(a, 5) - 0.251289 * pow(a, 6)) * pow(rd, 4)) * pow(rv, 2) +
				(0.0838975 + 0.198853 * a - 0.0994766 * pow(a, 2) - 0.822409 * pow(a, 3) -
					0.586006 * pow(a, 4) + 1.17402 * pow(a, 5) + 1.17378 * pow(a, 6) +
					(-0.182047 - 0.245637 * a + 0.130396 * pow(a, 2) + 2.85439 * pow(a, 3) +
						2.44277 * pow(a, 4) - 5.36406 * pow(a, 5) - 5.22806 * pow(a, 6)) * rd +
					(-0.227493 - 0.470604 * a - 0.670555 * pow(a, 2) - 1.66997 * pow(a, 3) +
						0.268677 * pow(a, 4) + 3.71822 * pow(a, 5) + 2.10678 * pow(a, 6)) * pow(rd, 2) +
					(-0.13635 - 0.193843 * a + 0.626076 * pow(a, 2) - 1.55192 * pow(a, 3) -
						3.22512 * pow(a, 4) + 3.03851 * pow(a, 5) + 4.01364 * pow(a, 6)) * pow(rd, 3) +
					(-0.0262554 - 0.0391291 * a + 0.312858 * pow(a, 2) - 0.122063 * pow(a, 3) -
						1.03112 * pow(a, 4) + 0.28978 * pow(a, 5) + 0.878604 * pow(a, 6)) * pow(rd, 4)) * pow(rv, 3))
				* (1. - rd);

	v = min(5., v);
//	cout << "pes v = " << v << endl;
//	getchar();
	return v;
}
