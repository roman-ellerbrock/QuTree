#include "NOCl.h"

namespace Operator {
// Library with new operators
	Tensorcd Applyqd(const LeafInterface& grid, const Tensorcd& Acoeffs, int exponent) {
		const Vectord& x = grid.getX();
		Tensorcd xAcoeffs(Acoeffs);

		// mu stuff
		const double mh = 1836.;
		const double mN = 14.0 * mh;
		const double mCl = 35.403 * mh;
		const double mO = 16.0 * mh;
//		const double muv = 1. / mN + 1. / mO;
//		const double mud = 1. / (mN + mO) + 1. / mCl;
		const double muv = 1./13615.5;
		const double mud = 1./29456.5;

//		cout << "muv: " << 1./muv << endl;
//		cout << "mud: " << 1./mud << endl;
//		cout << "muN: " << 1./mN << endl;
//		getchar();
		// equilibrium geometry
		const double alpha = 1.5;
		const double rde = 4.315;
		for (int i = 0; i < Acoeffs.shape().lastBefore(); i++) {
			double factor = 1 - exp(-alpha * (x(i) * sqrt(mud) - rde));
			factor = pow(factor, exponent);
			for (int n = 0; n < Acoeffs.shape().lastDimension(); n++) {
				xAcoeffs(i, n) *= factor;
			}
		}
		return xAcoeffs;
	}

// dissociation (R)
	void qd1(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqd(grid, Acoeffs, 1);
	}

	void qd2(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqd(grid, Acoeffs, 2);
	}

	void qd3(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqd(grid, Acoeffs, 3);
	}

	void qd4(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqd(grid, Acoeffs, 4);
	}

	void qdminus(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		Tensorcd xAcoeffs(Acoeffs.shape());
		qd1(grid, xAcoeffs, Acoeffs);
		hA = (Acoeffs - xAcoeffs);
	}

// vibration (r)
	Tensorcd Applyqv(const LeafInterface& grid, const Tensorcd& Acoeffs, int exponent) {
		const Vectord& x = grid.getX();
		Tensorcd xAcoeffs(Acoeffs);
		// mu stuff
		const double mh = 1836.;
		const double mN = 14.0 * mh;
		const double mCl = 35.403 * mh;
		const double mO = 16.0 * mh;
//		const double muv = 1. / mN + 1. / mO;
//		const double mud = 1. / (mN + mO) + 1. / mCl;
		const double muv = 1./13615.5;
		const double mud = 1./29456.5;
//		cout << "muv: " << 1./muv << endl;
//		cout << "mud: " << 1./mud << endl;
//		cout << "muN: " << mN << endl;
//		getchar();
		// equilibriu geometrys
		const double rve = 2.136;
		for (int i = 0; i < Acoeffs.shape().lastBefore(); i++) {
			double factor = pow(x(i) * sqrt(muv) - rve, exponent);
			for (int n = 0; n < Acoeffs.shape().lastDimension(); n++) {
				xAcoeffs(i, n) *= factor;
			}
		}
		return xAcoeffs;
	}

	void qv1(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqv(grid, Acoeffs, 1);
	}

	void qv2(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqv(grid, Acoeffs, 2);
	}

	void qv3(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqv(grid, Acoeffs, 3);
	}

	void qv4(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {
		hA = Applyqv(grid, Acoeffs, 4);
	}

// Parts for the moments of inertia
	void recrsq(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeffs) {

		const Vectord& x = grid.getX();
		Tensorcd xAcoeffs(Acoeffs);
		for (int n = 0; n < Acoeffs.shape().lastDimension(); n++) {
			for (int i = 0; i < Acoeffs.shape().lastBefore(); i++) {
				xAcoeffs(i, n) /= x(i) * x(i);
			}
		}
		hA = xAcoeffs;
	}

	void ApplyW00(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			// actual applying the operator here
			double factor = 0;
			factor += 0.0384816;
			factor += 0.0247875 * q;
			factor += 0.0270933 * q * q;
			factor += 0.00126791 * q * q * q;
			factor += 0.00541285 * q * q * q * q;
			factor += 0.0313629 * q * q * q * q * q;
			factor += 0.0172449 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW01(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.00834237;
			factor += 0.00398713 * q;
			factor += 0.00783319 * q * q;
			factor += 0.0294887 * q * q * q;
			factor += -0.0154387 * q * q * q * q;
			factor += -0.0621984 * q * q * q * q * q;
			factor += -0.0337951 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW02(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.00161625;
			factor += -0.00015633 * q;
			factor += -0.0189982 * q * q;
			factor += -0.00753297 * q * q * q;
			factor += 0.00383665 * q * q * q * q;
			factor += -0.00758225 * q * q * q * q * q;
			factor += -0.00904493 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW03(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.0010101;
			factor += 0.000619148 * q;
			factor += -0.0149812 * q * q;
			factor += -0.0199722 * q * q * q;
			factor += 0.00873013 * q * q * q * q;
			factor += 0.0376118 * q * q * q * q * q;
			factor += 0.0221523 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW04(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.0003689;
			factor += 0.000164037 * q;
			factor += -0.00331809 * q * q;
			factor += -0.00567787 * q * q * q;
			factor += 0.00268662 * q * q * q * q;
			factor += 0.0134483 * q * q * q * q * q;
			factor += 0.0084585 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW10(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.0558666;
			factor += -0.0276576 * q;
			factor += 0.0934932 * q * q;
			factor += -0.0295638 * q * q * q;
			factor += -0.15436 * q * q * q * q;
			factor += 0.0796119 * q * q * q * q * q;
			factor += 0.135121 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW11(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.0582169;
			factor += 0.0384404 * q;
			factor += 0.078114 * q * q;
			factor += 0.185556 * q * q * q;
			factor += -0.0641656 * q * q * q * q;
			factor += -0.175976 * q * q * q * q * q;
			factor += -0.0104994 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW12(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.052285;
			factor += 0.0472724 * q;
			factor += -0.216008 * q * q;
			factor += -0.147775 * q * q * q;
			factor += 0.349283 * q * q * q * q;
			factor += 0.28458 * q * q * q * q * q;
			factor += 0.00384449 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW13(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.0212609;
			factor += 0.0290597 * q;
			factor += -0.109124 * q * q;
			factor += 0.0310445 * q * q * q;
			factor += 0.262513 * q * q * q * q;
			factor += -0.250653 * q * q * q * q * q;
			factor += -0.369466 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW14(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.00334178;
			factor += 0.0039061 * q;
			factor += -0.0110452 * q * q;
			factor += 0.0582029 * q * q * q;
			factor += 0.0679524 * q * q * q * q;
			factor += -0.16459 * q * q * q * q * q;
			factor += -0.165337 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW20(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.163186;
			factor += -0.180535 * q;
			factor += 0.04692 * q * q;
			factor += 0.471673 * q * q * q;
			factor += 0.403267 * q * q * q * q;
			factor += -0.718071 * q * q * q * q * q;
			factor += -0.761199 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW21(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.0290674;
			factor += -0.0136172 * q;
			factor += -0.108952 * q * q;
			factor += -1.68269 * q * q * q;
			factor += -1.2673 * q * q * q * q;
			factor += 3.17648 * q * q * q * q * q;
			factor += 2.92793 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW22(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.121228;
			factor += 0.202308 * q;
			factor += 0.483613 * q * q;
			factor += 1.29095 * q * q * q;
			factor += -0.174483 * q * q * q * q;
			factor += -2.4605 * q * q * q * q * q;
			factor += -1.36597 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW23(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.107233;
			factor += 0.115213 * q;
			factor += -0.366102 * q * q;
			factor += 0.812662 * q * q * q;
			factor += 1.76038 * q * q * q * q;
			factor += -1.19665 * q * q * q * q * q;
			factor += -1.77392 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW24(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.0232767;
			factor += 0.0304932 * q;
			factor += -0.19455 * q * q;
			factor += -0.037517 * q * q * q;
			factor += 0.539365 * q * q * q * q;
			factor += 0.120203 * q * q * q * q * q;
			factor += -0.251289 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW30(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += 0.0838975;
			factor += 0.198853 * q;
			factor += -0.0994766 * q * q;
			factor += -0.822409 * q * q * q;
			factor += -0.586006 * q * q * q * q;
			factor += 1.17402 * q * q * q * q * q;
			factor += 1.17378 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW31(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.182047;
			factor += -0.245637 * q;
			factor += 0.130396 * q * q;
			factor += 2.85439 * q * q * q;
			factor += 2.44277 * q * q * q * q;
			factor += -5.36406 * q * q * q * q * q;
			factor += -5.22806 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW32(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.227493;
			factor += -0.470604 * q;
			factor += -0.670555 * q * q;
			factor += -1.66997 * q * q * q;
			factor += 0.268677 * q * q * q * q;
			factor += 3.71822 * q * q * q * q * q;
			factor += 2.10678 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW33(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.13635;
			factor += -0.193843 * q;
			factor += 0.626076 * q * q;
			factor += -1.55192 * q * q * q;
			factor += -3.22512 * q * q * q * q;
			factor += 3.03851 * q * q * q * q * q;
			factor += 4.01364 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyW34(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		const double beta = 1.1;
		const double rthetae = -0.60737584;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double q = exp(-beta * x(i)) - exp(-beta * rthetae);

			double factor = 0;
			// actual applying the operator here
			factor += -0.0262554;
			factor += -0.0391291 * q;
			factor += 0.312858 * q * q;
			factor += -0.122063 * q * q * q;
			factor += -1.03112 * q * q * q * q;
			factor += 0.28978 * q * q * q * q * q;
			factor += 0.878604 * q * q * q * q * q * q;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

	void ApplyNOpotential(const LeafInterface& grid, Tensorcd& hA, const Tensorcd& Acoeff) {
		const Vectord& x = grid.getX();
		// mu stuff
		const double mh = 1836.;
		const double mN = 14.0 * mh;
		const double mCl = 35.403 * mh;
		const double mO = 16.0 * mh;
//		const double muv = 1. / mN + 1. / mO;
//		const double mud = 1. / (mN + mO) + 1. / mCl;
		const double muv = 1./13615.5;
		const double mud = 1./29456.5;
		// equilibriu geometrys
		const double rve = 2.136;

		Tensorcd hAcoeff(Acoeff.shape());

		for (int i = 0; i < x.dim(); i++) {
			const double qv = x(i) * sqrt(muv) - rve;

			double factor = 0;
			// actual applying the operator here
			factor += 0.6816 * qv * qv;
			factor += -0.9123 * qv * qv * qv;
			factor += 0.4115 * qv * qv * qv * qv;
			for (int n = 0; n < Acoeff.shape().lastDimension(); n++) {
				hAcoeff(i, n) = factor * Acoeff(i, n);
			}
		}
		hA = hAcoeff;
	}

/*	Vectord NOCl::getW(vector<double> coeffs, const LeafInterface& grid)
	{
		const Vectord& x = grid.getX();
		const int dim = x.dim();
		Vectord W(dim);

		const double beta = 1.1;
		const double rthetae = -0.60737584;

		// new coordinate
		Vectord rtheta(dim);
		for (int i = 0; i < dim; i++) {
			rtheta(i) = exp(-beta*x(i)) - exp(-beta*rthetae);
		}

		// calculate contribution to the overall vector
		for (int n = 0; n < 6; n++) {
			for (int i = 0; i < dim; i++) {
				W(i) += coeffs[i] * pow(rtheta(i), n);
			}
		}
		return W;
	}*/

	SOPcd NOCl_V() {
		LeafFuncd qv1a = qv1;
		LeafFuncd qv2a = qv2;
		LeafFuncd qv3a = qv3;

		LeafFuncd qdminusa = qdminus;

		LeafFuncd qd1a = qd1;
		LeafFuncd qd2a = qd2;
		LeafFuncd qd3a = qd3;
		LeafFuncd qd4a = qd4;

		LeafFuncd ApplyW00a = ApplyW00;
		LeafFuncd ApplyW01a = ApplyW01;
		LeafFuncd ApplyW02a = ApplyW02;
		LeafFuncd ApplyW03a = ApplyW03;
		LeafFuncd ApplyW04a = ApplyW04;

		LeafFuncd ApplyW10a = ApplyW10;
		LeafFuncd ApplyW11a = ApplyW11;
		LeafFuncd ApplyW12a = ApplyW12;
		LeafFuncd ApplyW13a = ApplyW13;
		LeafFuncd ApplyW14a = ApplyW14;

		LeafFuncd ApplyW20a = ApplyW20;
		LeafFuncd ApplyW21a = ApplyW21;
		LeafFuncd ApplyW22a = ApplyW22;
		LeafFuncd ApplyW23a = ApplyW23;
		LeafFuncd ApplyW24a = ApplyW24;

		LeafFuncd ApplyW30a = ApplyW30;
		LeafFuncd ApplyW31a = ApplyW31;
		LeafFuncd ApplyW32a = ApplyW32;
		LeafFuncd ApplyW33a = ApplyW33;
		LeafFuncd ApplyW34a = ApplyW34;

		LeafFuncd ApplyNOpotentiala = ApplyNOpotential;

		SOPcd V;
		for (int i = 0; i <= 3; i++) {
			for (int j = 0; j <= 4; j++) {
				MLOcd M;

				// (1-q_d)
				M.push_back(qdminus, 1);

				// r
				if (i == 1) { M.push_back(qv1a, 0); }
				if (i == 2) { M.push_back(qv2a, 0); }
				if (i == 3) { M.push_back(qv3a, 0); }

				// R
				if (j == 1) { M.push_back(qd1a, 1); }
				if (j == 2) { M.push_back(qd2a, 1); }
				if (j == 3) { M.push_back(qd3a, 1); }
				if (j == 4) { M.push_back(qd4a, 1); }

				// theta
				if (i == 0 && j == 0) { M.push_back(ApplyW00a, 2); }
				if (i == 0 && j == 1) { M.push_back(ApplyW01a, 2); }
				if (i == 0 && j == 2) { M.push_back(ApplyW02a, 2); }
				if (i == 0 && j == 3) { M.push_back(ApplyW03a, 2); }
				if (i == 0 && j == 4) { M.push_back(ApplyW04a, 2); }

				if (i == 1 && j == 0) { M.push_back(ApplyW10a, 2); }
				if (i == 1 && j == 1) { M.push_back(ApplyW11a, 2); }
				if (i == 1 && j == 2) { M.push_back(ApplyW12a, 2); }
				if (i == 1 && j == 3) { M.push_back(ApplyW13a, 2); }
				if (i == 1 && j == 4) { M.push_back(ApplyW14a, 2); }

				if (i == 2 && j == 0) { M.push_back(ApplyW20a, 2); }
				if (i == 2 && j == 1) { M.push_back(ApplyW21a, 2); }
				if (i == 2 && j == 2) { M.push_back(ApplyW22a, 2); }
				if (i == 2 && j == 3) { M.push_back(ApplyW23a, 2); }
				if (i == 2 && j == 4) { M.push_back(ApplyW24a, 2); }

				if (i == 3 && j == 0) { M.push_back(ApplyW30a, 2); }
				if (i == 3 && j == 1) { M.push_back(ApplyW31a, 2); }
				if (i == 3 && j == 2) { M.push_back(ApplyW32a, 2); }
				if (i == 3 && j == 3) { M.push_back(ApplyW33a, 2); }
				if (i == 3 && j == 4) { M.push_back(ApplyW34a, 2); }

				// save operator
				V.push_back(M, 1.);
			}
		}

		{
			MLOcd M;
			M.push_back(ApplyNOpotentiala, 0);
			V.push_back(M, 1.);
		}
		return V;
	}

	SOPcd NOCl_KEO() {
		auto x = &LeafInterface::applyX;
		auto kin = &LeafInterface::applyKin;
		auto p = &LeafInterface::applyP;
		LeafFuncd recrsqa = recrsq;

		SOPcd H;
		// kinetic energy
		// r
		{
			MLOcd M;
			M.push_back(kin, 0);
			H.push_back(M, 1.);
		}
		// R
		{
			MLOcd M;
			M.push_back(kin, 1);
			H.push_back(M, 1.);
		}

		// Theta
		{
			MLOcd M;
			M.push_back(recrsqa, 0);
			M.push_back(kin, 2);
			H.push_back(M, 1.);
		}
		{
			MLOcd M;
			M.push_back(recrsqa, 1);
			M.push_back(kin, 2);
			H.push_back(M, 1.);
		}
		return H;
	}

	SOPcd NOCl_H(bool potential) {
		if (potential) {
			return NOCl_KEO() + NOCl_V();
		} else {
			return NOCl_KEO();
		}
	}
}