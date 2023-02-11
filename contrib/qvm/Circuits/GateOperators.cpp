//
// Created by Roman Ellerbrock on 2019-04-19.
//
#include "GateOperators.h"
#include "Core/Matrix_Implementation.h"
#include "TreeOperators/SumOfProductsOperator_Implementation.h"

namespace Circuits {

	void Hadamard(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = 1. / sqrt(2.) * (phi(0, n) + phi(1, n));
			xPhi(1, n) = 1. / sqrt(2.) * (phi(0, n) - phi(1, n));
		}
	}

	MLOcd Hadamard_mpo(const Register& r) {
		MLOcd M;
		for (size_t i = r.front(); i < r.end(); ++i) {
			M.push_back(Hadamard, i);
		}
		return M;
	}

	void X(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = phi(1, n);
			xPhi(1, n) = phi(0, n);
		}
	}

/*	Matrixcd X(bool adj) {
		Matrixcd x(2, 2);
		x(0, 0) = 0;
		x(1, 0) = 1;
		x(0, 1) = 1;
		x(1, 1) = 0;
		return x;
	}
 */

	LeafMatrixcd Xalpha(double alpha, bool adj) {
		Matrixcd x(2, 2);
		alpha *= QM::pi / 2.;
		x(0, 0) = cos(alpha);
		x(1, 0) = -QM::im * sin(alpha);
		x(0, 1) = -QM::im * sin(alpha);
		x(1, 1) = cos(alpha);
		x *= exp(QM::im * alpha);
		return LeafMatrixcd(x, adj);
	}

	LeafMatrixcd Zalpha(double alpha, bool adj) {
		Matrixcd z(2, 2);
		z(0, 0) = 1.;
		z(1, 1) = exp(QM::im * QM::pi * alpha);
		return LeafMatrixcd(z, adj);
	}

	void Y(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = -QM::im * phi(1, n);
			xPhi(1, n) = QM::im * phi(0, n);
		}
	}

	void Z(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = -phi(1, n);
		}
	}

	void S(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = QM::im * phi(1, n);
		}
	}

	void Sdagger(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = -QM::im * phi(1, n);
		}
	}

	void set0(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = 0.;
		}
	}

	void set1(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = 0;
			xPhi(1, n) = phi(1, n);
		}
	}

	void identity(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		for (size_t i = 0; i < phi.shape().totalDimension(); ++i) {
			xPhi(i) = phi(i);
		}
	}

	void SetBlack(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorShape& tdim = xPhi.shape();
//		double epsilon = 0.95;
		double epsilon = 0.10;
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			xPhi(0, n) = epsilon * phi(0, n) + sqrt(1. - pow(epsilon, 2)) * phi(1, n);
			xPhi(1, n) = sqrt(1. - pow(epsilon, 2)) * phi(0, n) + epsilon * phi(0, n);
		}
	}

	Matrixcd T() {
		Matrixcd T(2, 2);
		T(0, 0) = 1.;
		T(1, 1) = (1. + QM::im) / (sqrt(2.));
		return T;
	}

	// https://arxiv.org/pdf/quant-ph/0511250.pdf
	Matrixcd sqrtX() {
		Matrixcd sX(2, 2);
		complex<double> coeff = 1. / 2.;
		sX(0, 0) = 1. + QM::im;
		sX(0, 1) = 1. - QM::im;
		sX(1, 0) = 1. - QM::im;
		sX(1, 1) = 1. + QM::im;
		sX *= coeff;
		return sX;
	}

	// https://arxiv.org/pdf/quant-ph/0511250.pdf
	Matrixcd sqrtY() {
		Matrixcd sY(2, 2);
		complex<double> coeff = 1. / 2.;
		sY(0, 0) = 1. + QM::im;
		sY(0, 1) = -1. - QM::im;
		sY(1, 0) = 1. + QM::im;
		sY(1, 1) = 1. + QM::im;
		sY *= coeff;
		return sY;
	}

	/// Define complex rotation operator, denoted by Rk, where k indicates angle=1/(2^k)
	Rk::Rk(size_t k, bool adjungate) {
		long double pw = pow(2, k);
		complex<long double> im(0., 1.);
		long double two_pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640 * 2.;
		complex<long double> phase(im * two_pi / pw);
		if (adjungate) { phase = conj(phase); }
		factor = exp(phase);
	}

	void Rk::apply(const LeafInterface& grid, Tensorcd& RPhi, const Tensorcd& Phi) const {
		const TensorShape& tdim = RPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			RPhi(0, n) = Phi(0, n);
			RPhi(1, n) = factor * Phi(1, n);
		}
	}

	/// Define complex rotation operator, denoted by Rx, where x indicates the angle
	Rx::Rx(long double x, bool adjungate) {
		complex<long double> im(0., 1.);
		long double two_pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640 * 2.;
		complex<long double> phase(im * two_pi * x);
		if (adjungate) { phase = conj(phase); }
		factor = exp(phase);
	}

	void Rx::apply(const LeafInterface& grid, Tensorcd& RPhi, const Tensorcd& Phi) const {
		const TensorShape& tdim = RPhi.shape();
		for (size_t n = 0; n < tdim.lastDimension(); ++n) {
			RPhi(0, n) = Phi(0, n);
			RPhi(1, n) = factor * Phi(1, n);
		}
	}

	Matrixcd sigma_x() {
		Matrixcd s(2, 2);
		s(0, 0) = 0.;
		s(1, 0) = 1.;
		s(0, 1) = 1.;
		s(1, 1) = 0.;
		return s;
	}

	Matrixcd sigma_y() {
		Matrixcd s(2, 2);
		s(0, 0) = 0.;
		s(1, 0) = QM::im;
		s(0, 1) = -QM::im;
		s(1, 1) = 0.;
		return s;
	}

	Matrixcd sigma_z() {
		Matrixcd s(2, 2);
		s(0, 0) = 1.;
		s(1, 0) = 0.;
		s(0, 1) = 0.;
		s(1, 1) = -1.;
		return s;
	}

	Matrixcd rX(double theta) {
		Matrixcd s(2, 2);
		s(0, 0) = cos(theta / 2.);
		s(1, 0) = -QM::im * sin(theta / 2.);
		s(0, 1) = -QM::im * sin(theta / 2.);
		s(1, 1) = cos(theta / 2.);
		return s;
	}

	Matrixcd rY(double theta) {
		Matrixcd s(2, 2);
		s(0, 0) = cos(theta / 2.);
		s(1, 0) = sin(theta / 2.);
		s(0, 1) = sin(theta / 2.);
		s(1, 1) = cos(theta / 2.);
		return s;
	}

	Matrixcd rZ(double theta) {
		Matrixcd s(2, 2);
		s(0, 0) = exp(-QM::im * theta / 2.);
		s(1, 0) = 0.;
		s(0, 1) = 0.;
		s(1, 1) = exp(QM::im * theta / 2.);
		return s;
	}

	Matrixcd U(double theta, double phi, double lambda) {
		Matrixcd s(2, 2);
		s(0, 0) = cos(theta / 2.);
		s(1, 0) = exp(QM::im * phi) * sin(theta / 2.);
		s(0, 1) = -exp(QM::im * lambda) * sin(theta / 2.);
		s(1, 1) = exp(QM::im * (phi + lambda)) * cos(theta / 2.);
		return s;
	}

    Matrixcd givens_rot(double theta) {
        Matrixcd s(2, 2);
        s(0, 0) = cos(theta);
        s(1, 0) = -sin(theta);
        s(0, 1) = sin(theta);
        s(1, 1) = cos(theta);
        return s;
    }

    Matrixcd sigmaPlus() {
        /// Ref. [1] Eq. (19)
        Matrixcd s(2, 2);
        s(1, 0) = 1.;
        return s;
//    return 0.5 * (sigmaX() - QM::im * sigmaY());
    }

    Matrixcd sigmaMinus() {
        /// Ref. [1] Eq. (19)
        Matrixcd s(2, 2);
        s(0, 1) = 1.;
        return s;
//    return 0.5 * (sigmaX() + QM::im * sigmaY());
    }


	SOPVectorcd ising(size_t q0, size_t q1, size_t steps) {
		double rz1 = 1. / (double) steps;
		double rz2 = rz1 / 4.;
		double rx = rz1;
		SOPVectorcd I;
		I.append(CNot(q0, q1));
		MLOcd M(LeafMatrixcd(rZ(rz1)), q1);
		I.append(M);
		I.append(CNot(q0, q1));
		auto Rs = rX(rx) * rZ(rz2);
		MLOcd rs;
		rs.push_back(Rs, q0);
		rs.push_back(Rs, q1);
		I.append(rs);
		return I;
	}

	/* Follow conventions in arXiv:2106.13839 */
	// returns SOP with six nonzero terms: 1, c, -s, s, c, 1
	SOPcd givens(double theta, size_t c, size_t t){ // c for control, not cosine
	    SOPcd S;
        {
            MLOcd M;
            M.push_back(set0,c);
            M.push_back(set0,t);
            S.push_back(M,1);
        }
        {
            MLOcd M;
            M.push_back(set1,c);
            M.push_back(set1,t);
            S.push_back(M,1);
        }
        {
            MLOcd M;
            M.push_back(set0,c);
            M.push_back(set1,t);
            S.push_back(M,cos(theta));
        }
        {
            MLOcd M;
            M.push_back(set1,c);
            M.push_back(set0,t);
            S.push_back(M,cos(theta));
        }
        {
            MLOcd M;
            M.push_back(sigmaMinus(),c);
            M.push_back(sigmaPlus(),t);
            S.push_back(M,-sin(theta));
        }
        {
            MLOcd M;
            M.push_back(sigmaPlus(),c);
            M.push_back(sigmaMinus(),t);
            S.push_back(M,sin(theta));
        }
        return S;

	}

	// CGivens is a controlled Givens rotation on a subspace of qubits q1 and q2
    SOPcd cGivens(double theta, size_t control, size_t q1, size_t q2) {
        return makeCGate(givens(theta, q1, q2), control);
    }

	LeafFunctioncd randomProject(mt19937& gen) {
		uniform_real_distribution<double> dist(0., 1.);
		if (dist(gen) > 0.5) {
			return LeafFunctioncd(set0);
		} else {
			return LeafFunctioncd(set1);
		}
	}

	SOPcd makeCGate(const SOPcd& in, size_t control) {
		using namespace Circuits;
		SOPcd out;
		{
			// Apply identity to 0 projection
			MLOcd M;
			M.push_back(set0, control);

			out.push_back(M, 1.0);
		}

		{
			// Apply given gate to 1 projection
			MLOcd M;
			M.push_back(set1, control);

			// Distribute through existing SOP
			SOPcd Cin = M * in;

			// Concatenate with Set0
			out = out + Cin;
		}

		return out;
	}

	/**
	 * LeafOperator h = X
	 * MultiLeafOperator (MLOcd) M = X(0) * Y(1) * ...
	 * SumOfProductOperator (SOPcd) S = 1.0 * M + 2.0 * N + ...
	 */

	SOPcd makeCCGate(const SOPcd& in, size_t control1, size_t control2) {
		using namespace Circuits;
		SOPcd out;

		{
			// Apply identity if first control is 0
			// Catches both second control cases
			MLOcd M;
			M.push_back(set0, control1);
			out.push_back(M, 1.0);
		}

		{
			// Apply identity if second control is 0
			MLOcd M;
			M.push_back(set1, control1);
			M.push_back(set0, control2);
			out.push_back(M, 1.0);
		}

		{
			// Apply given gate if both controls are 1
			MLOcd M;
			M.push_back(set1, control1);
			M.push_back(set1, control2);

			SOPcd Cin = M * in;

			out = out + Cin;
		}

		return out;
	}

//	CNot = Set0(c)*1(t) + Set1(c)*X(t)
	SOPcd CNot_expanded(size_t control, size_t target) {
		MLOcd M1(set0, control);

		MLOcd M2(set1, control);
		M2.push_back(X, target);

		SOPcd S;
		S.push_back(M1, 1.);
		S.push_back(M2, 1.);
		return S;
	}

	SOPcd CNot(size_t control, size_t target) {
		return makeCGate(X, control, target);
	}

	SOPcd CNot(const Register& r1, const Register& r2) {
		return CNot(r1.front(), r2.front());
	}

	SOPcd CZ(size_t control, size_t target) {
		return makeCGate(Z, control, target);
	}

	SOPcd CZ(const Register& r1, const Register& r2) {
		return CZ(r1.front(), r2.front());
	}

	SOPcd Toffoli(size_t control1, size_t control2, size_t target) {
		return makeCCGate(X, control1, control2, target);
	}

	SOPcd Toffoli(const Register& c1, const Register& c2, const Register& t) {
		return Toffoli(c1.front(), c2.front(), t.front());
	}

	SOPcd Uk(size_t control, size_t target, size_t k, bool conj_phase) {
		shared_ptr<LeafOperatorcd> R = make_shared<Rk>(k, conj_phase);
		return makeCGate(R, control, target);
	}

	SOPcd Uk(const Register& c, const Register& t, size_t k, bool conj_phase) {
		return Uk(c.front(), t.front(), k, conj_phase);
	}

	SOPcd Ux(size_t c, size_t t, double x, bool conj_phase) {
		shared_ptr<LeafOperatorcd> R = make_shared<Rx>(x, conj_phase);
		return makeCGate(R, c, t);
	}

	SOPcd Ux(const Register& c, const Register& t, double x, bool conj_phase) {
		return Ux(c.front(), t.front(), x, conj_phase);
	}

	SOPcd CANzeroSOP(size_t c, size_t t, double alphax, double alphay, bool adj) {
		MLOcd Zs(Z, c);
		Zs.push_back(Z, t);
		MLOcd Vs(sqrtX(), c);
		Vs.push_back(sqrtX(), t);
		MLOcd XZ(Xalpha(alphax, adj), c);
		XZ.push_back(Zalpha(alphay, adj), t);

		SOPcd canzero(Zs);
		canzero = Vs * canzero;
		canzero = CNot(c, t) * canzero;
		canzero = XZ * canzero;
		canzero = CNot(c, t) * canzero;
		canzero = Zs * canzero;
		canzero = Vs * canzero;
		return canzero;
	}

	SOPVectorcd CANzero(size_t c, size_t t, double alphax, double alphay, bool adj) {
		SOPVectorcd canzero;
		MLOcd Zs(Z, c);
		Zs.push_back(Z, t);
		MLOcd Vs(sqrtX(), c);
		Vs.push_back(sqrtX(), t);

		MLOcd XZ(Xalpha(alphax, adj), c);
		XZ.push_back(Zalpha(alphay, adj), t);

		canzero.append(Zs);
		canzero.append(Vs);
		canzero.append(CNot(c, t));
		canzero.append(XZ);
		canzero.append(CNot(c, t));
		canzero.append(Zs);
		canzero.append(Vs);
		return canzero;
	}

	SOPcd iSWAPSOP(const size_t c, const size_t t, double theta, bool adj) {
		auto circ = CANzeroSOP(c, t, 0.25, 0.25, adj);
		circ = Ux(c, t, theta / 2., adj) * circ;
		return circ;
	}

	SOPVectorcd iSWAP(const size_t c, const size_t t, double theta, bool adj) {
		auto circ = CANzero(c, t, 0.25, 0.25, adj);
		circ.append(Ux(c, t, theta / 2., adj));
		return circ;
	}

	SOPVectorcd distribSingleControl(const SOPVectorcd& stack, size_t control) {
		SOPVectorcd Cstack;

		for (int i = 0; i < stack.size(); ++i) {
			SOPcd Cgate = makeCGate(stack[i], control);
			Cstack.append(Cgate);
		}

		return Cstack;
	}

	SOPVectorcd distribDoubleControl(const SOPVectorcd& stack, size_t control1, size_t control2) {
		SOPVectorcd CCstack;

		for (int i = 0; i < stack.size(); ++i) {
			SOPcd CCgate = makeCCGate(stack[i], control1, control2);
			CCstack.append(CCgate);
		}

		return CCstack;
	}

	SOPVectorcd Carry(size_t carry1, size_t a, size_t b, size_t carry2) {
		/*! \brief CARRY operation is a basic circuit in quantum computing.
		 *
		 * The circuit looks like this (* = control, x = flip)
		 * c1 -----------*---
		 * a  ---*---*-------
		 * b  ---*---x---*---
		 * c2 ---x-------x---
		 */
		SOPVectorcd stack;
		stack.emplace_back(Toffoli(a, b, carry2));
		stack.emplace_back(CNot(a, b));
		stack.emplace_back(Toffoli(carry1, b, carry2));
		return stack;
	}

	SOPVectorcd Carry(const Register& c1, const Register& a, const Register& b, const Register& c2) {
		return Carry(c1.front(), a.front(), b.front(), c2.front());
	}

	SOPVectorcd RCarry(size_t carry1, size_t a, size_t b, size_t carry2) {
		/*! \brief reverse CARRY operation is a basic circuit in quantum computing.
		 *
		 * The circuit is the Carry operation in reverse order
		 */
		SOPVectorcd stack = Carry(carry1, a, b, carry2);
		std::reverse(stack.begin(), stack.end());
		return stack;
	}

	SOPVectorcd RCarry(const Register& c1, const Register& a, const Register& b, const Register& c2) {
		return RCarry(c1.front(), a.front(), b.front(), c2.front());
	}

	SOPVectorcd SUM(size_t a, size_t b, size_t sum) {
		SOPVectorcd sops;
		sops.emplace_back(CNot(a, sum));
		sops.emplace_back(CNot(b, sum));
		return sops;
	}

	SOPVectorcd SUM(const Register& a, const Register& b, const Register& sum) {
		return SUM(a.front(), b.front(), sum.front());
	}

	SOPVectorcd CNotChain(const Register& reg) {
		SOPVectorcd stack;
		for (size_t i = reg.front(); i < reg.end() - 1; ++i) {
			stack.push_back(CNot(i, i + 1));
		}
		return stack;
	}

	SOPVectorcd HadamardChain(const Register& reg) {
		MLOcd Hs;
		for (size_t k = reg.front(); k < reg.end(); ++k) {
			Hs.push_back(Hadamard, k);
		}
		SOPVectorcd stack;
		stack.append(SOPcd(Hs, 1.));
		return stack;
	}

	SOPVectorcd distribute(const Register& reg, function<void(const LeafInterface&, Tensorcd&, const Tensorcd&)> f) {
		MLOcd M;
		for (size_t k = reg.front(); k < reg.end(); ++k) {
			M.push_back(f, k);
		}
		SOPVectorcd stack;
		stack.append(SOPcd(M, 1.));
		return stack;
	}

	SOPVectorcd CNotBenchmark(const Register& reg) {
		SOPVectorcd stack;
		{
			MLOcd M(Hadamard, reg.front());
			SOPcd sop(M, 1.);
			stack.push_back(sop);
		}
		stack.append(CNotChain(reg));
		return stack;
	}

	SOPVectorcd Adder_circuit(const Tree& tree, size_t n_bit,
		size_t mode_a, size_t mode_b, size_t mode_carry) {
		/*! \brief This circuit adds two numbers
		 *
		 *
		 */

		size_t f = tree.nLeaves();
		SOPVectorcd stack;
		// Perform the carry part of the summation
		for (size_t n = 0; n < n_bit; ++n) {
			size_t c1 = mode_carry + n;
			size_t a = mode_a + n;
			size_t b = mode_b + n;
			size_t c2 = mode_carry + n + 1;
			assert (c1 < f && a < f && b < f && c2 < f);
			cout << "Carry(" << c1 << ", " << a << ", " << b << ", " << c2 << ")\n";
			auto carry = Carry(c1, a, b, c2);
			stack.append(carry);
//			stack.insert(stack.end(), carry.begin(), carry.end());
		}

		// A single CNot operation
		{
			size_t a = mode_a + n_bit - 1;
			size_t b = mode_b + n_bit - 1;
			cout << "CNot(" << a << ", " << b << ")\n";
			stack.emplace_back(CNot(a, b));
		}

		// Perform
		for (int n = n_bit - 1; n >= 0; --n) {
			size_t c1 = mode_carry + n;
			size_t a = mode_a + n;
			size_t b = mode_b + n;
			size_t c2 = mode_carry + n + 1;
			assert (c1 < f && a < f && b < f && c2 < f);
			if (n != (n_bit - 1)) {
				cout << "RCarry(" << c1 << ", " << a << ", " << b << ", " << c2 << ")\n";
				auto rcarry = RCarry(c1, a, b, c2);
				stack.insert(stack.end(), rcarry.begin(), rcarry.end());
			}
			cout << "SUM(" << c1 << ", " << a << ", " << b << ")\n";
			auto sum = SUM(c1, a, b);
			stack.insert(stack.end(), sum.begin(), sum.end());
		}
		return stack;
	}

	SOPVectorcd RandomSuperposition(const Tree& tree, size_t start, size_t end, size_t n) {
		assert(end <= tree.nLeaves());
		assert(n <= (end - start));

		srand(time(nullptr));
		int seed = rand();
		mt19937 gen(seed);
		vector<size_t> arr;
		for (int i = start; i < end; ++i)
			arr.push_back(i);
		std::shuffle(arr.begin(), arr.end(), gen);
		MLOcd Hs;
		for (size_t k = 0; k < n; ++k) {
			cout << "mode=" << arr[k] << endl;
			assert(arr[k] < tree.nLeaves());
			Hs.push_back(Hadamard, arr[k]);
		}
		SOPVectorcd stack;
		stack.append(SOPcd(Hs, 1.));
		return stack;
	}

	size_t Idx(size_t ix, size_t iy, size_t n_row) {
		return iy * n_row + ix;
	}

	Matrix<bool> true_mat(size_t n_row, size_t n_col) {
		Matrix<bool> m(n_row, n_col);
		for (size_t row = 0; row < n_row; ++row) {
			for (size_t col = 0; col < n_col; ++col) {
				m(row, col) = true;
			}
		}
		return m;
	}

	SOPVectorcd CNotGridRow(const Tree& tree,
		size_t n_row, size_t n_col, size_t offset, size_t shift) {
		Matrix<bool> m = true_mat(n_row, n_col);

		SOPVectorcd stack;
		for (size_t row = 0; row < n_row; ++row) {
			size_t offset_row = 2 * ((row + offset) % 2);
			for (size_t col = shift; col < n_col - offset_row - shift; col += 4) {
				size_t eff_col = col + offset_row;
				size_t idx1 = Idx(row, eff_col, n_row);
				size_t idx2 = Idx(row, eff_col + 1, n_row);
//				assert((idx1 < n_row) && (idx2 < n_col));
				cout << idx1 << " " << idx2 << endl;
				cout << "( " << row << ", " << eff_col << ") - ( " << row << ", " << eff_col + 1 << ")\n";
				m(row, eff_col) = false;
				m(row, eff_col + 1) = false;
				idx1 = ConvertIdx2D(idx1, n_row, n_col);
				idx2 = ConvertIdx2D(idx2, n_row, n_col);
				stack.emplace_back(CNot(idx1, idx2));
			}
		}
		m.print();
		return stack;
	}

	SOPVectorcd CNotGridCol(const Tree& tree, size_t n_row, size_t n_col, size_t offset, size_t shift) {
		Matrix<bool> m = true_mat(n_row, n_col);
		m.print();

		SOPVectorcd stack;
		for (size_t col = 0; col < n_col; col++) {
			size_t offset_col = 2 * ((col + offset) % 2);
			for (size_t row = shift; row < n_row - offset_col - shift; row += 4) {
				size_t eff_row = row + offset_col;
				size_t idx1 = Idx(eff_row, col, n_row);
				size_t idx2 = Idx(eff_row + 1, col, n_row);
				m(eff_row, col) = false;
				m(eff_row + 1, col) = false;
				cout << "( " << eff_row << ", " << col << ") - ( " << eff_row + 1 << ", " << col << ")\n";
//				assert((idx1 < n_row) && (idx2 < n_col));
				idx1 = ConvertIdx2D(idx1, n_row, n_col);
				idx2 = ConvertIdx2D(idx2, n_row, n_col);
				stack.emplace_back(CNot(idx1, idx2));
			}
		}
		m.print();
		return stack;
	}

	SOPVectorcd CNotGridPattern(const Tree& tree, size_t n_row, size_t n_col, size_t which) {
		switch (which) {
			case 0:return CNotGridRow(tree, n_row, n_col, 0, 0);
			case 1:return CNotGridRow(tree, n_row, n_col, 0, 1);
			case 2:return CNotGridRow(tree, n_row, n_col, 1, 0);
			case 3:return CNotGridRow(tree, n_row, n_col, 1, 1);
			case 4:return CNotGridCol(tree, n_row, n_col, 0, 0);
			case 5:return CNotGridCol(tree, n_row, n_col, 0, 1);
			case 6:return CNotGridCol(tree, n_row, n_col, 1, 0);
			case 7:return CNotGridCol(tree, n_row, n_col, 1, 1);
			default:cerr << "Grid pattern are defined for patterns 0,..,7\n";
				exit(1);
		}
	}

	SOPVectorcd CNotGridSwipe(const Tree& tree, size_t n_row, size_t n_col) {
		SOPVectorcd stack;
		for (size_t n = 0; n < 8; ++n) {
			cout << "n = " << n << endl;
			SOPVectorcd GridOp = CNotGridPattern(tree, n_row, n_col, n);
			stack.append(GridOp);
		}
		return stack;
	}

	SOPVectorcd GoogleIterations(const Tree& tree, size_t n_row, size_t n_col, size_t n_swipes) {
		SOPVectorcd stack;
		{
			MLOcd M;
			SOPcd sop;
			for (size_t n = 0; n < tree.nLeaves(); n += 2) {
				M.push_back(Hadamard, n);
			}
			sop.push_back(M, 1.0);
			stack.push_back(sop);
		}
		for (size_t n = 0; n < n_swipes; ++n) {
			auto SwipeStack = CNotGridSwipe(tree, n_row, n_col);
			stack.insert(stack.end(), SwipeStack.begin(), SwipeStack.end());
		}
		return stack;
	}

	size_t ConvertIdx2D(size_t idx, size_t nrow, size_t ncol) {
		/** \brief This routine maps a linear index to a suitable index in a binary quadratic lattice
		 *
		 * Mapping:
		 * 0  1  | 2  3      0  1  | 4  5
		 * ----- | ----      ----- | ----
		 * 4  5  | 6  7      2  3  | 6  7
		 * ------------- ==> -------------
		 * 8  9  | 10 11     8  9  | 12 13
		 * ----- | ----      ----- | ----
		 * 12 13 | 14 15     10 11 | 14 15
		 *
		 * This routine will map a index of the map left to the index of the map on the right hand side.
		 * n is the number of indices per coordinate (4 in this example). nrow & ncol must be powers of 2.
		 * */
//        assert(nrow == ncol);
		size_t col = idx % ncol;
		size_t row = idx / ncol;
		size_t new_idx = 0;
		size_t nprod = nrow * ncol;
		bool horizontal;
		while (nprod > 0) {
			(nrow >= ncol) ? horizontal = true : horizontal = false;
			nprod /= 2;
			if (horizontal) {
				nrow /= 2;
//                cout << "row" << row << ">?" << nrow << " +" << nprod << " =" << new_idx << endl;
				if (row >= nrow) {
					row -= nrow;
					new_idx += nprod;
				}
			} else {
				ncol /= 2;
//                cout << "col" << col << ">?" << ncol << " +" << nprod << " =" << new_idx << endl;
				if (col >= ncol) {
					col -= ncol;
					new_idx += nprod;
				}
			}
		}
		return new_idx;
	}

	SOPVectorcd distribute(const LeafFunPaircd& h, const Register& reg, bool adjoint) {
		MLOcd M;
		for (size_t i = 0; i < reg.size(); ++i) {
			size_t idx = reg.front() + i;
			M.push_back(h, idx, adjoint);
		}
		return SOPVectorcd(SOPcd(M));
	}

	SOPVectorcd distribute(const LeafFuncd& h, const Register& reg, bool adjoint) {
		return distribute({h, h}, reg, adjoint);
	}

	SOPVectorcd distibute(const Matrixcd& h, const Register& reg, bool adjoint) {
		MLOcd M;
		for (size_t i = 0; i < reg.size(); ++i) {
			size_t idx = reg.front() + i;
			M.push_back(h, idx, adjoint);
		}
		return SOPVectorcd(SOPcd(M));
	}

	SOPVectorcd SWAP(const size_t c, const size_t t) {
		SOPVectorcd circ({CNot(c, t)});
		circ.emplace_back(CNot(t, c));
		circ.emplace_back(CNot(c, t));
		return circ;
	}
}
