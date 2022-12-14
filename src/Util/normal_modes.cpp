//
// Created by Roman Ellerbrock on 4/28/22.
//

#include "Util/normal_modes.h"
#include "TreeOperators/CoordinateTransformation.h"

double Vmassweighted(Vectord q, const Potential& V, const Matrixd& U) {
	q = U * q;
	return V.evaluate(q, 0);
}

/// in cartesian coordinates
Vectord grad(const Vectord& x, const Potential& V, const Vectord& h, const Matrixd& U) {
	Vectord G(x.dim());
	for (size_t k = 0; k < G.dim(); ++k) {

		Vectord xp = x;
		xp(k) += h[k];
		double ve = Vmassweighted(xp, V, U);

		Vectord xm = x;
		xm(k) -= h[k];
		double va = Vmassweighted(xm, V, U);

		G(k) = (ve - va) / (2. * h[k]);
	}
	return G;
}

/// in internal coordinates
Vectord grad(const Vectord& x, const PotentialOperator& V, double delta) {
	Vectord G(x.dim());
	for (size_t k = 0; k < G.dim(); ++k) {

		Vectord xe = x;
		xe(k) += delta;
		double ve = V.evaluate(xe, 0);

		Vectord xa = x;
		xa(k) -= delta;
		double va = V.evaluate(xa, 0);
		G(k) = (ve - va) / (2. * delta);
	}
	return G;
}

Vectord getq0(const Tree& tree) {
	Vectord q(tree.nLeaves());
	for (size_t l = 0; l < tree.nLeaves(); ++l) {
		const Leaf& leaf = tree.getLeaf(l);
		q(leaf.mode()) = leaf.par().wfr0();
	}
	return q;
}

Matrixd hessian(const Vectord& X, const Potential& V,
	const Vectord& h, const Matrixd& U) {

	vector<Vectord> grad_plus;
	vector<Vectord> grad_minus;
	for (size_t n = 0; n < X.dim(); ++n) {
		{
			Vectord Xp(X);
			Xp(n) += h[n];
			grad_plus.emplace_back(grad(Xp, V, h, U));
		}
		{
			Vectord Xm(X);
			Xm(n) -= h[n];
			grad_minus.emplace_back(grad(Xm, V, h, U));
		}
	}

	Matrixd H(X.dim(), X.dim());
	for (size_t j = 0; j < X.dim(); ++j) {
		for (size_t i = 0; i < X.dim(); ++i) {
			H(i, j) =
				(grad_plus[j](i) - grad_minus[j](i)) / (4 * h[j])
					+ (grad_plus[i](j) - grad_minus[i](j)) / (4 * h[i]);
		}
	}
	for (size_t i = 0; i < X.dim(); ++i) {
		double v = Vmassweighted(X, V, U);
		Vectord xp(X);
		xp(i) += h[i];
		double vp = Vmassweighted(xp, V, U);
		Vectord xpp(X);
		xpp(i) += 2 * h[i];
		double vpp = Vmassweighted(xpp, V, U);

		Vectord xm(X);
		xm(i) -= h[i];
		double vm = Vmassweighted(xm, V, U);
		Vectord xmm(X);
		xmm(i) -= 2 * h[i];
		double vmm = Vmassweighted(xmm, V, U);

		H(i, i) = (-vpp + 16 * vp - 30 * v + 16 * vm - vmm) / (12 * h[i] * h[i]);
	}
	return H;
}

void normalModeOutput(const Matrixd& Hessian, const Vectord& m, const Vectord& X0) {

	auto diag = diagonalize(Hessian);
	auto ev = diag.second;
	auto U = diag.first;

	const double ha2cm = 219474.6313705;
	cout << "==========================================\n";
	cout << "Normal Mode frequencies (cm-1):\n";
	cout << "==========================================\n";
	for (size_t i = 0; i < ev.dim(); ++i) {
		if (ev(i) > 0) {
			cout << i << "\t-\t" << sqrt(ev(i)) * ha2cm << " | " << sqrt(ev(i)) << endl;
		} else {
			cout << i << "\t-\t-" << sqrt(abs(ev(i))) * ha2cm << " | " << sqrt(abs(ev(i))) << endl;
		}
	}

	ofstream os("U.dat");
	ofstream os_x("x0.dat");
	os << U.dim1() << " " << U.dim2() << endl;
	for (size_t i = 0; i < U.dim1(); ++i) {
		for (size_t j = 0; j < U.dim2(); ++j) {
			os << U(j, i) / m(j) << " ";
		}
		os << endl;
	}
	for (size_t i = 0; i < U.dim1(); ++i) {
		os_x << X0(i) << " ";
	}
	double E = 0;
	for (size_t i = 6; i < U.dim1(); ++i) {
		E += 0.5 * sqrt(abs(ev(i)));
	}
	cout << "Harmonic ground state energy: " << E * ha2cm << " cm-1" << endl;
}

void find_minimum(const Hamiltonian& H, const Tree& tree) {
	const PotentialOperator& V = H.V_;
	Vectord q = getq0(tree);
/*	for (size_t i = 0; i < q.dim(); ++i) {
		q(i) += 0.01;
	}

	size_t n_iter = 500;
	double hlen = 1e-2;
	double delta = 1e-3;

	for (size_t iter = 0; iter < n_iter; iter++) {
		cout << "iter = " << iter << ", v = " << V.evaluate(q, 0) << endl;
		auto G = grad(q, V, delta);
		q -= G * hlen;
	}*/

	const Potential& v = *V.v();
	Vectord X = V.Q_->transform(q);

	cout << "Final energy: " << v.evaluate(X, 0) << endl;

	size_t natom = 4;
	double mh = 1836.153;

	Vectord mass_atom(natom * 3);
	mass_atom(0) = 12.00;
	mass_atom(1) = 1.00;
	mass_atom(2) = 1.00;
	mass_atom(3) = 1.00;
	Vectord h(natom * 3);
	Vectord m(natom * 3);
	double delta = 1e-3;
	Matrixd U(h.dim(), h.dim());
	for (size_t i = 0; i < 3 * natom; ++i) {
		m(i) = sqrt(mass_atom(i / 3) * mh);
		h(i) = delta * m(i);
		U(i, i) = 1. / m(i);
	}

	Vectord X0 = X;
	for (size_t i = 0; i < 3 * natom; ++i) {
		X(i) *= m(i);
	}

	cout << "Hessian at energy: " << Vmassweighted(X, v, U) << endl;

	auto Hessian = hessian(X, v, h, U);

	normalModeOutput(Hessian, m, X0);
}

