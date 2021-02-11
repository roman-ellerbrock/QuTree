#include "Util/WeightedSimultaneousDiagonalization.h"

namespace WeightedSimultaneousDiagonalization {

	void Calculate(vector<Matrixcd>& Xs, vector<Matrixcd> XXs,
		Matrixcd& W, Matrixcd& trafo, double eps) {
		// Checks
		for (const Matrixcd& x : Xs) {
			assert(W.dim1() == x.dim1());
		}
		for (const Matrixcd& x : XXs) {
			assert(W.dim1() == x.dim1());
		}
		assert(Xs.size() > 0);

		// Set Initial variables
		int iter = 0;
		int maxiter = 25;
		bool converged = false;

		// Initialize the Diagonalization by rotating to the diagonal
		// representation of one of the matrices in Xs. This
		// avoids stationary points during the optimization process.
		trafo.zero();
		for (size_t i = 0; i < trafo.dim1(); i++)
			trafo(i, i) = 1.;

		// Initial transformation of the matrices
		vector<Matrixcd> Xs_plain(Xs);
		WeightMatrices(Xs, W);
		WeightMatrices(XXs, W);

		// Iterate Jacobirotations until a converged result is reached
		// Measure off-diagonal norm
		double delta = MeasureWeightedOffDiagonality(Xs, Xs_plain, W, trafo);
//		cout << "Start : " << delta << endl;
		double delta_last = 2 * delta;
		while (!converged && iter < maxiter) {
			// Rotation circle over all elements
			WeightedJacobiRotations(Xs, XXs, W, trafo);

			// Measure off-diagonal norm
			delta_last = delta;
			delta = MeasureWeightedOffDiagonality(Xs, Xs_plain, W, trafo);
			double change = abs(delta_last - delta);
//			cout << iter << " " << delta << " " << change << endl;

			// Check convergence
			if (delta < eps || change < eps) { converged = true; }
			iter++;
		}

		if (!converged) {
//			cout << "D_WSD: " << MeasureWeightedOffDiagonality(Xs, Xs_plain, W, trafo) << endl;
		}
	}

	double MeasureWeightedDiagonality(
		const vector<Matrixcd>& A, const Matrixcd& W) {
		// Measure of diagonality in WSD
		double loc = 0;
		for (const Matrixcd& X : A) {
			// loc += - Xwii**2/Wii
			for (size_t i = 0; i < X.dim1(); i++) {
				double wdiag = real(W(i, i));
				double xdiag = real(X(i, i));
				loc += xdiag * xdiag / wdiag;
			}
		}
		return loc;
	}

	void WeightedJacobiRotations(
		vector<Matrixcd>& Xs, vector<Matrixcd>& XXs, Matrixcd& W, Matrixcd& trafo) {
		// Angles for Givens-Rotation
		complex<double> c, s = 0;

		// Swipe over the matrix-dimensions and perform jacobi-rotations
		for (size_t i = 0; i < W.dim1() - 1; i++) {
			for (size_t j = i + 1; j < W.dim1(); j++) {
				// Calculate Angles c and s for the elements i and j
				CalculateWeightedAngles(c, s, i, j, Xs, XXs, W);
//			cout << "c, s = " << c << " " << s << endl;

				assert(abs(1. - abs(c) * abs(c) - abs(s) * abs(s)) < 1E-10);

				// Perform the Givens-Rotation with angles c and s on the Xw-Matrices
				RotateMatrices(Xs, c, s, i, j);

				// Perform the Givens-Rotation with angles c and s on the Xw²-Matrices
				RotateMatrices(XXs, c, s, i, j);

				// The same for the Weight matrix
				GivensRotation(W, c, s, i, j);

				// Rotate the Transformation-Matrix
				GivensTrafoRotation(trafo, c, s, i, j);
			}
		}
	}

	int CalculateWeightedAngles(
		complex<double>& c, complex<double>& s,
		size_t i, size_t j, const vector<Matrixcd>& Xs,
		const vector<Matrixcd>& XXs, const Matrixcd& W) {
		// Rational function optimization (RFO) to optimize angles in WSD

		// Control parameters for RFO-Algorithm
		int iter = 0;
		int maxiter = 30;
		// Convergence factor and starting distortion of angles
		double eps = 1E-3;
		// arbitrary starting value for rel. diff.
		double diff_rel = 2 * eps * eps;
		// Set "a" for RFO
		double amin = 1E-12;
		double amax = 1;
		double a = amax;
		double delta = 1E-5;

		// Starting Angles for Search
		Vectord snew(2), cnew(2);
		Vectord svec(2), cvec(2);
		svec(0) = eps;
		svec(1) = -eps;
		cvec(0) = sqrt(1 - 2 * eps * eps);
		cvec(1) = 0;
		c = InterpretComplex(cvec);
		s = InterpretComplex(svec);

		// Calculate localization measure for the starting angles
		double loc_alt = WeightedJacobiLoc(Xs, XXs, W, i, j, c, s);
		double loc_new = 0;

		while (iter < maxiter && diff_rel > eps * eps && a > amin) {
			// Calculate derivatives of weightes measure
			Matrixd preHessian(2, 2);
			Vectord grad(2);
			WeightedJacobiDerivatives(grad, preHessian, s, Xs, XXs, W, i, j, delta);

			double diff = -1;
			// Iteration for Rational Function Optimization
			while (a > amin && diff < 0 && diff_rel > pow(eps, 2)) {
				double z = 1;
				Matrixd trafo(3, 3);
				// Decrease a, until z <= 0.5
				while (a > amin && z > 0.5) {
					// Build RFO-Hessian
					Matrixd H = RFO_BuildHessian(preHessian, grad, a);
					// Diagonalize and adjust eigenvectors
					Vectord ev(3);
					H.rDiag(trafo, ev);
					Vectord deltaS(2);
					deltaS(0) = a * trafo(0, 2) / trafo(2, 2);
					deltaS(1) = a * trafo(1, 2) / trafo(2, 2);

					// Calculate new angles
					snew(0) = real(s) + deltaS(0);
					snew(1) = imag(s) + deltaS(1);
					// z = |s_new|�
					complex<double> z_cd = InterpretComplex(snew);
					z = pow(abs(z_cd), 2);
					if (z > 0.5) { a /= 2; }
				}
				cnew(0) = sqrt(abs(1 - z));
				cnew(1) = 0;

				// Calculate Weighted loc for new angles
				complex<double> cnow = InterpretComplex(cnew);
				complex<double> snow = InterpretComplex(snew);
				loc_new = WeightedJacobiLoc(Xs, XXs, W, i, j, cnow, snow);

				// Check if the measure is decreased
				double trafoabs = pow(trafo(0, 2), 2) + pow(trafo(1, 2), 2);
				double sabs = pow(abs(s), 2) + eps * eps;
				diff_rel = trafoabs / sabs;
				diff = loc_new - loc_alt;

				// Adjust a-factor
				if (diff < 0) {
					a /= 2;
				} else {
					if (a < amax) { a *= 2; }
				}
			}

			// save angles, if they improve the measure
			if (a > amin && diff > 0) {
				c = InterpretComplex(cnew);
				s = InterpretComplex(snew);
				loc_alt = loc_new;
			}
			iter++;
			// iterate until diff_rel < eps*eps
		}

		// Check if Angle-optimization worked (0) or not (1)
		if (iter > maxiter || a > amin) {
			return 1;
		} else {
			return 0;
		}
	}

	/*
	double WeightedJacobiLoc(
			const vector<Matrixcd>& Xs, const vector<Matrixcd>& XXs,
			const Matrixcd& W, size_t p_, size_t q,
			complex<double> c, complex<double> s) {
		// Weight matrix contribution
		// First rotation W*J^H
		constexpr double lambda = 1e-4;
//		constexpr double lambda = 0.;
		Vectord W_new = RotatedDiagonals(W, p_, q, c, s);

		// Change of Xw-Matrix diagonals
		Vectord delta(2), correct(2);
		Vectord x_old(2), x_new(2), W_old(2);
		W_old(0) = real(W(p_, p_));
		W_old(1) = real(W(q, q));
		for(size_t i = 0; i < Xs.size(); ++i) {
			const Matrixcd& Xw = Xs[i];
			const Matrixcd& XXw = XXs[i];
			// Calculate the reduced measure for first order (x_)
			delta(0) = pow(real(Xw(p_, p_)), 2) / W_old(0);
			delta(1) = pow(real(Xw(q, q)), 2) / W_old(1);
			// Calculate the correction for the second order (x_²)
			correct(0) = lambda * (-2 * real(XXw(p_, p_)) + delta(0)) / W_old(0);
			correct(1) = lambda * (-2 * real(XXw(q, q)) + delta(1)) / W_old(1);
			x_old(0) += delta(0) * (1 + correct(0));
			x_old(1) += delta(1) * (1 + correct(1));
		}

		for(size_t i = 0; i < Xs.size(); ++i) {
			const Matrixcd& Xw = Xs[i];
			const Matrixcd& XXw = XXs[i];
			Vectord Xw_rot = RotatedDiagonals(Xw, p_, q, c, s);
			Vectord XXw_rot = RotatedDiagonals(XXw, p_, q, c, s);
			// Calculate the reduced measure for first order (x_)
			delta(0) = pow(Xw_rot(0), 2) / (W_new(0));
			delta(1) = pow(Xw_rot(1), 2) / (W_new(1));
			// Calculate the correction for the second order (x_²)
			correct(0) = lambda * (-2 * XXw_rot(0) + delta(0)) / W_new(0);
			correct(1) = lambda * (-2 * XXw_rot(1) + delta(1)) / W_new(1);
			x_new(0) += delta(0) * (1 + correct(0));
			x_new(1) += delta(1) * (1 + correct(1));
		}

		// Resulting change in the localization measure
		double diff = 0;
		diff += x_new(0) - x_old(0);
		diff += x_new(1) - x_old(1);

		return diff;
	}
	 */

	double WeightedJacobiLoc(
		const vector<Matrixcd>& Xs, const vector<Matrixcd>& XXs,
		const Matrixcd& W, size_t p, size_t q,
		complex<double> c, complex<double> s) {
		// Weight matrix contribution
		// First rotation W*J^H
		auto W_new = RotatedDiagonals(W, p, q, c, s);

		// Change of Xw-Matrix diagonals
//		Vectord x_old(2), x_new(2);
		double xold0 = 0, xold1 = 0, xnew0 = 0, xnew1 = 0;
		for (const Matrixcd& Xw : Xs) {
			xold0 += pow(real(Xw(p, p)), 2);
			xold1 += pow(real(Xw(q, q)), 2);
		}

		for (const Matrixcd& Xw : Xs) {
			auto Xw_rot = RotatedDiagonals(Xw, p, q, c, s);
			xnew0 += pow(Xw_rot.first, 2);
			xnew1 += pow(Xw_rot.second, 2);
		}

		// Resulting change in the localization measure
		double diff = 0;
		diff += xnew0 / W_new.first - xold0 / real(W(p, p));
		diff += xnew1 / W_new.second - xold1 / real(W(q, q));

		return diff;
	}

	void WeightedJacobiDerivatives(Vectord& grad,
		Matrixd& preHessian, complex<double> s_in,
		const vector<Matrixcd>& Xs, const vector<Matrixcd>& XXs,
		const Matrixcd& W, size_t p, size_t q, double delta) {
		assert(delta > 0);
		assert(grad.dim() == 2);
		assert(preHessian.dim1() == 2);
		assert(preHessian.dim2() == 2);
		Vectord s(2), x(2);
		s(0) = real(s_in);
		s(1) = imag(s_in);

		// Build a displacement-matrix y
		Matrixd y(3, 3);
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				x(0) = s(0) + ((int) i - 1) * delta;
				x(1) = s(1) + ((int) j - 1) * delta;
				double z = x(0) * x(0) + x(1) * x(1);
				complex<double> cnow(sqrt(abs(1 - z)), 0);
				complex<double> snow = InterpretComplex(x);
				y(i, j) = WeightedJacobiLoc(Xs, XXs, W, p, q, cnow, snow);
			}
		}

		// Calculate the gradient from the displacement matrix
		grad(0) = (y(2, 1) - y(0, 1)) / (2 * delta);
		grad(1) = (y(1, 2) - y(1, 0)) / (2 * delta);

		// Calculate the hessian from the displacement matrix
		preHessian(0, 0) = (y(2, 1) + y(0, 1) - 2 * y(1, 1)) / (delta * delta);
		preHessian(1, 1) = (y(1, 2) + y(1, 0) - 2 * y(1, 1)) / (delta * delta);
		preHessian(0, 1) = (y(2, 2) + y(0, 0) - y(2, 0) - y(0, 2)) / (4 * delta * delta);
		preHessian(1, 0) = preHessian(0, 1);
	}

	double MeasureWeightedOffDiagonality(
		const vector<Matrixcd>& Xws, const vector<Matrixcd>& Xs,
		const Matrixcd& W, const Matrixcd& trafo) {
		// Measure of diagonality in WSD
		double loc = 0;
		for (size_t k = 0; k < Xws.size(); k++) {
			const Matrixcd& Xw = Xws[k];
			const Matrixcd& X = Xs[k];
			Matrixcd X_diag(Xw.dim1(), Xw.dim2());

			for (size_t i = 0; i < Xw.dim1(); i++) {
				X_diag(i, i) = real(Xw(i, i) / W(i, i));
			}

			// Xd = trafo_^A * Xd * trafo_
			Matrixcd X_trafo = unitarySimilarityTrafo(X, trafo);

			// X-Xtrafo
			X_diag = X_diag - X_trafo;
			X_diag = X_diag * X_diag;

			// Weight with density matrix
			X_diag = W * X_diag;

			loc += real(X_diag.trace());
		}
		return loc;
	}

	vector<Vectord> GetDiagonals(const vector<Matrixcd>& Xws, const Matrixcd& W) {
		vector<Vectord> x_evs;
		for (const Matrixcd& x : Xws) {
			Vectord x_ev(x.dim1());
			for (size_t i = 0; i < x.dim1(); i++) {
				x_ev(i) = real(x(i, i) / W(i, i));
			}
			x_evs.emplace_back(x_ev);
		}

		return x_evs;
	}

	pair<Matrixcd, vector<Vectord>>
	Calculate(vector<Matrixcd>& Xs, Matrixcd& W, double eps) {
		auto trafo = identityMatrix<complex<double>>(W.dim1());
		vector<Matrixcd> XXs;

		Calculate(Xs, XXs, W, trafo, eps);

		vector<Vectord> x_evs = GetDiagonals(Xs, W);
		return {trafo, x_evs};
	}
}
