//
// Created by Roman Ellerbrock on 7/27/21.
//

#include "JordanWigner.h"
#include "Util/QMConstants.h"

namespace JordanWigner {
	Matrixcd sigmaX() {
		Matrixcd x(2, 2);
		x(1, 0) = 1.;
		x(0, 1) = 1.;
		return x;
	}

	Matrixcd sigmaY() {
		Matrixcd y(2, 2);
		y(1, 0) = QM::im;
		y(0, 1) = -QM::im;
		return y;
	}

	Matrixcd sigmaZ() {
		Matrixcd z(2, 2);
		z(0, 0) = 1.;
		z(1, 1) = -1.;
		return z;
	}

	Matrixcd sigmaPlus() {
		/// Ref. [1] Eq. (19)
		return 0.5 * (sigmaX() - QM::im * sigmaY());
	}

	Matrixcd sigmaMinus() {
		/// Ref. [1] Eq. (19)
		return 0.5 * (sigmaX() + QM::im * sigmaY());
	}

	MLOcd creation(size_t p) {
		/// Ref. [1], Eq. (17), M = (\prod_(i<p) sigmaZ_i) * sigmaMinus_p
		MLOcd M(sigmaPlus(), p);
		for (size_t i = 0; i < p; ++i) {
			M.push_back(sigmaZ(), i);
		}
		return M;
	}

	MLOcd annihilation(size_t p) {
		/// Ref. [1], Eq. (18), M = (\prod_(i<p) sigmaZ_i) * sigmaMinus_p
		MLOcd M(sigmaMinus(), p);
		for (size_t i = 0; i < p; ++i) {
			M.push_back(sigmaZ(), i);
		}
		return M;
	}

	MLOcd twoIndexOperator(size_t p, size_t q, double eps) {
		size_t max = p + 1;
		if (q >= max) max = q + 1;

		MLOcd A; /// A = (a_p^+) (a_q)
		cout << p << " " << q << " | ";
		for (size_t i = 0; i < max; ++i) {
			Matrixcd op = identityMatrixcd(2);

			if (i < q) {
				op = sigmaZ() * op;
				cout << "z^q_" << i << " ";
			} else if (i == q) {
				op = sigmaMinus() * op;
				cout << "a^q_" << i << " ";
			}

			if (i < p) {
				op = sigmaZ() * op;
				cout << "z^p_" << i << " ";
			} else if (i == p) {
				op = sigmaPlus() * op;
				cout << "a^p+_" << i << " ";
			}

			if ((residual(op, identityMatrixcd(2)) > eps) || (i == (max - 1))) {
				A.push_back(op, i);
				cout << "! ";
			}
		}
		cout << "| size = " << A.size() << endl;
		return A;
	}

	MLOcd fourIndexOperator(size_t p, size_t q, size_t r, size_t s, double eps) {
		size_t max = p + 1;
		if (q >= max) max = q + 1;
		if (r >= max) max = r + 1;
		if (s >= max) max = s + 1;

		cout << p << " " << q << " " << r << " " << s << " | ";
		MLOcd A; /// A = (a_p^+) (a_q^+) (a_r) (a_s)
		/// Rational: swipe over each
		for (size_t i = 0; i < max; ++i) {
			Matrixcd op = identityMatrixcd(2);
			if (i < s) {
				op = sigmaZ() * op;
				cout << "z^s_" << i << " ";
			} else if (i == s) {
				op = sigmaMinus() * op;
				cout << "a^s_" << i << " ";
			}

			if (i < r) {
				op = sigmaZ() * op;
				cout << "z^r_" << i << " ";
			} else if (i == r) {
				op = sigmaMinus() * op;
				cout << "a^r_" << i << " ";
			}

			if (i < q) {
				op = sigmaZ() * op;
				cout << "z^q_" << i << " ";
			} else if (i == q) {
				op = sigmaPlus() * op;
				cout << "a^q+_" << i << " ";
			}

			if (i < p) {
				op = sigmaZ() * op;
				cout << "z^p_" << i << " ";
			} else if (i == p) {
				op = sigmaPlus() * op;
				cout << "a^p+_" << i << " ";
			}

			if (op.frobeniusNorm() < eps) {
			//	cout << "SMALL!\n";
			//	continue;
			}
			if ((residual(op, identityMatrixcd(2)) > eps) || (i == (max - 1))) {
				A.push_back(op, i);
				cout << "! ";
			}
		}
		cout << " | size = " << A.size() << endl;
		return A;
	}

	SOPcd electronicHamiltonian(const TwoIndex& Hpq, const FourIndex& Hpqrs) {
		double eps = 1e-10;
		SOPcd H;

		cout << "==== Hpq =================" << endl;
		for (const auto& hpq : Hpq) {
			size_t p = get<0>(hpq);
			size_t q = get<1>(hpq);
			double h = get<2>(hpq);
			if (abs(h) < eps) { continue; }

			auto M = twoIndexOperator(p, q, eps);
			if (M.size() == 0) { continue; }
			if (M.size() == 0) { cerr << "M empty in Hpq!\n"; exit(1); }
			H.push_back(M, h);
		}
		cout << "Electronic Hamiltonian size (Hpq): " << H.size() << endl;

		cout << "==== Hpqrs =================" << endl;
		for (const auto& term : Hpqrs) {
			size_t p = get<0>(term);
			size_t q = get<1>(term);
			size_t r = get<2>(term);
			size_t s = get<3>(term);
			double h = get<4>(term);
			if (abs(h) < eps) { continue; }

//			if (p == q || r == s) { continue; } /// always evaluates to zero due to a_p+ a_p+=0 or a_r a_r=0
			auto M = fourIndexOperator(p, q, r, s, eps);
			if (M.size() == 0) { continue; }
			if (M.size() == 0) { cerr << "M empty in Hpqrs!\n"; exit(1); }
			H.push_back(M, 0.5 * h);
		}
		cout << "Electronic Hamiltonian size (total): " << H.size() << endl;
		return H;
	}
}

