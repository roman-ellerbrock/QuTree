//
// Created by Roman Ellerbrock on 8/27/21.
//

#include "Util/GateOperators.h"
#include "Util/QMConstants.h"


Tensorcd sigma_x() {
	Tensorcd s({2, 2});
	s(0, 0) = 0.;
	s(1, 0) = 1.;
	s(0, 1) = 1.;
	s(1, 1) = 0.;
	return s;
}

Tensorcd sigma_y() {
	Tensorcd s({2, 2});
	s(0, 0) = 0.;
	s(1, 0) = QM::im;
	s(0, 1) = -QM::im;
	s(1, 1) = 0.;
	return s;
}

Tensorcd sigma_z() {
	Tensorcd s({2, 2});
	s(0, 0) = 1.;
	s(1, 0) = 0.;
	s(0, 1) = 0.;
	s(1, 1) = -1.;
	return s;
}

Tensorcd hadamard() {
	Tensorcd h({2, 2});
	h(0, 0) = 1. / sqrt(2.);
	h(1, 0) = 1. / sqrt(2.);
	h(0, 1) = 1. / sqrt(2.);
	h(1, 1) = -1. / sqrt(2.);
	return h;
}

Tensorcd element(size_t i, size_t j, size_t dim) {
	Tensorcd P({dim, dim});
	if (i >= dim || j >= dim) {
		cerr << "Error: element out of range in element-operator\n";
		exit(1);
	}
	P(i, j) = 1.;
	return P;
}

Tensorcd project(size_t i, size_t dim) {
	return element(i, i, dim);
}

SOPcd controlled(const SOPcd& in, size_t control) {
	SOPcd out;
	{
		// Apply identity to 0 projection
		POcd M;
		M.push_back(project(0, 2), control);

		out.push_back(M, 1.0);
	}

	{
		// Apply given gate to 1 projection
		POcd M;
		M.push_back(project(1, 2), control);

		// Distribute through existing SOP
		SOPcd Cin = M * in;

		// Concatenate with Set0
		out = out + Cin;
	}

	return out;
}

SOPcd controlled(const POcd& M, size_t c) {
	SOPcd S(M, 1.);
	return controlled(S, c);
}

SOPcd controlled(const Tensorcd& L, size_t c, size_t t) {
	POcd M(L, t);
	return controlled(M, c);
}

SOPcd CNot(size_t c, size_t t) {
	SOPcd cnot;
	{
		POcd M(project(0, 2), c);
		cnot.push_back(M, 1.);
	}
	{
		POcd M(project(1, 2), c);
		M.push_back(sigma_x(), t);
		cnot.push_back(M, 1.);
	}
	return cnot;
//	return controlled(sigma_x(), c, t);
}
/// CNot Chain
// (P(0;1) + P(1;1)*X(2)) * (P(0;0) + P(1;0)*X(1))
// (P(0;0)*P(0;1) +  P(0;0)*P(1;1)*X(2) +
// P(1;0)*X(1)*P(0;1) + P(1;0)*X(1)*P(1;1)*X(2))
// (second line only) = P(1;0)*P(1;1) + P(1;0)*P(0;1)*X(2)
// = P(0;0)*P(0;1) +  P(0;0)*P(1;1)*X(2) + P(1;0)*P(1;1) + P(1;0)*P(0;1)*X(2)
// = (P(0;0)*P(0;1) + P(1;0)P(1;1)) + (P(0;0)*P(1;1) + P(1;0) * P(0;1)) * X(2)
// P(0; 1)


SOPcd CZ(size_t c, size_t t) {
	return controlled(sigma_z(), c, t);
}

/// Define complex rotation operator, denoted by Rk, where k indicates angle=1/(2^k)
Tensorcd rk(size_t k, bool adj) {
	long double pw = pow(2, k);
	complex<long double> im(0., 1.);
	long double two_pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640 * 2.;
	complex<long double> phase(im * two_pi / pw);
	if (adj) { phase = -phase; }
	Tensorcd rmat({2, 2});
	rmat(0, 0) = 1.;
	rmat(1, 1) = exp(phase);
	return rmat;
}

SOPcd Uk(size_t c, size_t t, size_t k, bool adj) {
	return controlled(rk(k, adj), c, t);
}

/*SOPVectorcd QFT(size_t mode_start, size_t n_bit, bool adjungate, size_t approx) {
	SOPVectorcd qft;
	size_t mode_end = mode_start + n_bit;
	for (size_t n = 0; n < n_bit; ++n) {
		size_t bit_a = mode_start + n;
		{
			POcd M(hadamard(), bit_a);
			SOPcd sop(M, 1.);
			qft.push_back(sop);
		}
		//			cout << "H(" << bit_a << ")\n";
		for (size_t j = n + 1; j < n_bit; ++j) {
			size_t angle_k = j + 1 - n;
			size_t bit_b = mode_start + j;
			if (approx == 0 || angle_k < approx) {
				//					cout << "U(" << bit_a << ", " << bit_b << ")[" << angle_k << "] " << endl;
				qft.emplace_back(Uk(bit_a, bit_b, angle_k, adjungate));
			}
		}
	}
	return qft;
}*/

