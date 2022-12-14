//
// Created by Roman Ellerbrock on 2019-04-19.
//
#include "GateOperators.h"

namespace GateOperators {

	void Hadamard(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = 1. / sqrt(2.) * (phi(0, n) + phi(1, n));
			xPhi(1, n) = 1. / sqrt(2.) * (phi(0, n) - phi(1, n));
		}
	}

	MPO Hadamard_mpo(const Register& r) {
		MPO M;
		for (size_t i = r.Begin(); i < r.End(); ++i) {
			M.push_back(Hadamard, i);
		}
		return M;
	}

	void X(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = phi(1, n);
			xPhi(1, n) = phi(0, n);
		}
	}

	void Y(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = -QM::im * phi(1, n);
			xPhi(1, n) = QM::im * phi(0, n);
		}
	}

	void Z(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = -phi(1, n);
		}
	}

	void S(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = QM::im * phi(1, n);
		}
	}

	void T(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = (1. + QM::im) / (sqrt(2.)) * phi(1, n);
		}
	}

	void Sdagger(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = -QM::im * phi(1, n);
		}
	}

	void Set0(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = phi(0, n);
			xPhi(1, n) = 0.;
		}
	}

	void Set1(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = 0;
			xPhi(1, n) = phi(1, n);
		}
	}

	void Identity(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		for (size_t i = 0; i < phi.Dim().getdimtot(); ++i) {
			xPhi(i) = phi(i);
		}
	}

	void SetBlack(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
//		double epsilon = 0.95;
		double epsilon = 0.10;
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = epsilon * phi(0, n) + sqrt(1. - pow(epsilon, 2)) * phi(1, n);
			xPhi(1, n) = sqrt(1. - pow(epsilon, 2)) * phi(0, n) + epsilon * phi(0, n);
		}
	}

	// https://arxiv.org/pdf/quant-ph/0511250.pdf
	void sqrtX(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		complex<double> coeff = 1. / (2. * QM::im);
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = coeff * (QM::im * phi(0, n) + 1. * phi(1, n));
			xPhi(1, n) = coeff * (1. * phi(0, n)     + QM::im * phi(1, n));
		}
	}

	// https://arxiv.org/pdf/quant-ph/0511250.pdf
	void sqrtY(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi) {
		const TensorDim& tdim = xPhi.Dim();
		complex<double> coeff = 1. / (2. * QM::im);
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			xPhi(0, n) = coeff * (QM::im * phi(0, n) - QM::im * phi(1, n));
			xPhi(1, n) = coeff * (QM::im * phi(0, n) + QM::im * phi(1, n));
		}
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

	void Rk::Apply(Tensorcd& RPhi, const PrimitiveBasis& grid, const Tensorcd& Phi) const {
		const TensorDim& tdim = RPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
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

	void Rx::Apply(Tensorcd& RPhi, const PrimitiveBasis& grid, const Tensorcd& Phi) const {
		const TensorDim& tdim = RPhi.Dim();
		for (size_t n = 0; n < tdim.getntensor(); ++n) {
			RPhi(0, n) = Phi(0, n);
			RPhi(1, n) = factor * Phi(1, n);
		}
	}

	SOP MakeCGate(const SOP& in, size_t control) {
		using namespace GateOperators;
		SOP out;
		{
			// Apply identity to 0 projection
			MPO M;
			M.push_back(Set0, control);

			out.push_back(M, 1.0);
		}

		{
			// Apply given gate to 1 projection
			MPO M;
			M.push_back(Set1, control);

			// Distribute through existing SOP
			SOP Cin = M * in;

			// Concatenate with Set0
			out = out + Cin;
		}

		return out;
	}

	SOP MakeCCGate(const SOP& in, size_t control1, size_t control2) {
		using namespace GateOperators;
		SOP out;

		{
			// Apply identity if first control is 0
			// Catches both second control cases
			MPO M;
			M.push_back(Set0, control1);
			out.push_back(M, 1.0);
		}

		{
			// Apply identity if second control is 0
			MPO M;
			M.push_back(Set1, control1);
			M.push_back(Set0, control2);
			out.push_back(M, 1.0);
		}

		{
			// Apply given gate if both controls are 1
			MPO M;
			M.push_back(Set1, control1);
			M.push_back(Set1, control2);

			SOP Cin = M * in;

			out = out + Cin;
		}

		return out;
	}

//	CNot = Set0(c)*1(t) + Set1(c)*X(t)
	SOP CNot_expanded(size_t control, size_t target) {
		MPO M1(Set0, control);

		MPO M2(Set1, control);
		M2.push_back(X, target);

		SOP S;
		S.push_back(M1, 1.);
		S.push_back(M2, 1.);
		return S;
	}

	SOP CNot(size_t control, size_t target) {
		return MakeCGate(X, control, target);
	}

	SOP CNot(const Register& r1, const Register& r2) {
		return CNot(r1.Begin(), r2.Begin());
	}

	SOP Toffoli(size_t control1, size_t control2, size_t target) {
		return MakeCCGate(X, control1, control2, target);
	}

	SOP Toffoli(const Register& c1, const Register& c2, const Register& t) {
		return Toffoli(c1.Begin(), c2.Begin(), t.Begin());
	}

	SOP Uk(size_t control, size_t target, size_t k, bool conj_phase) {
		shared_ptr<BottomLayerSPO> R = make_shared<Rk>(k, conj_phase);
		return MakeCGate(R, control, target);
	}

	SOP Uk(const Register& c, const Register& t, size_t k, bool conj_phase) {
		return Uk(c.Begin(), t.Begin(), k, conj_phase);
	}

	SOPVector DistribSingleControl(const SOPVector& stack, size_t control) {
		SOPVector Cstack;

		for (int i = 0; i < stack.size(); ++i) {
			SOP Cgate = MakeCGate(stack[i], control);
			Cstack.append(Cgate);
		}

		return Cstack;
	}

	SOPVector DistribDoubleControl(const SOPVector& stack, size_t control1, size_t control2) {
		SOPVector CCstack;

		for (int i = 0; i < stack.size(); ++i) {
			SOP CCgate = MakeCCGate(stack[i], control1, control2);
			CCstack.append(CCgate);
		}

		return CCstack;
	}

	SOPVector Carry(size_t carry1, size_t a, size_t b, size_t carry2) {
		/*! \brief CARRY operation is a basic circuit in quantum computing.
		 *
		 * The circuit looks like this (* = control, x = flip)
		 * c1 -----------*---
		 * a  ---*---*-------
		 * b  ---*---x---*---
		 * c2 ---x-------x---
		 */
		SOPVector stack;
		stack.emplace_back(Toffoli(a, b, carry2));
		stack.emplace_back(CNot(a, b));
		stack.emplace_back(Toffoli(carry1, b, carry2));
		return stack;
	}

	SOPVector Carry(const Register& c1, const Register& a, const Register& b, const Register& c2) {
		return Carry(c1.Begin(), a.Begin(), b.Begin(), c2.Begin());
	}

	SOPVector RCarry(size_t carry1, size_t a, size_t b, size_t carry2) {
		/*! \brief reverse CARRY operation is a basic circuit in quantum computing.
		 *
		 * The circuit is the Carry operation in reverse order
		 */
		SOPVector stack = Carry(carry1, a, b, carry2);
		std::reverse(stack.begin(), stack.end());
		return stack;
	}

	SOPVector RCarry(const Register& c1, const Register& a, const Register& b, const Register& c2) {
		return RCarry(c1.Begin(), a.Begin(), b.Begin(), c2.Begin());
	}

	SOPVector SUM(size_t a, size_t b, size_t sum) {
		SOPVector sops;
		sops.emplace_back(CNot(a, sum));
		sops.emplace_back(CNot(b, sum));
		return sops;
	}

	SOPVector SUM(const Register& a, const Register& b, const Register& sum) {
		return SUM(a.Begin(), b.Begin(), sum.Begin());
	}

	SOPVector CNotChain(const Register& reg) {
		SOPVector stack;
		for (size_t i = reg.Begin(); i < reg.End() - 1; ++i) {
			stack.push_back(CNot(i, i + 1));
		}
		return stack;
	}

	SOPVector HadamardChain(const Register& reg) {
		MPO Hs;
		for (size_t k = reg.Begin(); k < reg.End(); ++k) {
			Hs.push_back(Hadamard, k);
		}
		SOPVector stack;
		stack.append(SOP(Hs, 1.));
		return stack;
	}

	SOPVector Distribute(const Register& reg, function<void(const PrimitiveBasis&, Tensorcd&, const Tensorcd&)> f) {
		MPO M;
		for (size_t k = reg.Begin(); k < reg.End(); ++k) {
			M.push_back(f, k);
		}
		SOPVector stack;
		stack.append(SOP(M, 1.));
		return stack;
	}

	SOPVector CNotBenchmark(const Register& reg) {
		SOPVector stack;
		{
			MultiParticleOperator M(Hadamard, reg.Begin());
			SOP sop(M, 1.);
			stack.push_back(sop);
		}
		stack.append(CNotChain(reg));
		return stack;
	}

	SOPVector Adder_circuit(const mctdhBasis& basis, size_t n_bit,
		size_t mode_a, size_t mode_b, size_t mode_carry) {
		/*! \brief This circuit adds two numbers
		 *
		 *
		 */

		size_t f = basis.nLeaves();
		SOPVector stack;
		// Perform the carry part of the summation
		for (size_t n = 0; n < n_bit; ++n) {
			size_t c1 = mode_carry + n;
			size_t a = mode_a + n;
			size_t b = mode_b + n;
			size_t c2 = mode_carry + n + 1;
			assert (c1 < f && a < f && b < f && c2 < f);
			cout << "Carry(" << c1 << ", " << a << ", " << b << ", " << c2 << ")\n";
			auto carry = Carry(c1, a, b, c2);
			stack.insert(stack.end(), carry.begin(), carry.end());
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

	SOPVector ControlledSuperposition(const Register& reg, size_t n) {
		MPO Hs;
		for (size_t i = reg.Begin(); i < reg.End(); ++i) {
			if ((i + 1) % n == 0) {
				cout << "mode=" << i << endl;
				Hs.push_back(Hadamard, i);
			}
		}
		SOPVector stack;
		stack.append(SOP(Hs, 1.));
		return stack;
	}

	SOPVector RandomSuperposition(const mctdhBasis& basis, size_t start, size_t end, size_t n) {
		assert(end <= basis.nLeaves());
		assert(n <= (end - start));

		srand(time(nullptr));
		int seed = rand();
		mt19937 gen(seed);
		vector<size_t> arr;
		for (int i = start; i < end; ++i)
			arr.push_back(i);
		std::shuffle(arr.begin(), arr.end(), gen);
		MPO Hs;
		for (size_t k = 0; k < n; ++k) {
			cout << "mode=" << arr[k] << endl;
			assert(arr[k] < basis.nLeaves());
			Hs.push_back(Hadamard, arr[k]);
		}
		SOPVector stack;
		stack.append(SOP(Hs, 1.));
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

	SOPVector CNotGridRow(const mctdhBasis& basis,
		size_t n_row, size_t n_col, size_t offset, size_t shift) {
		Matrix<bool> m = true_mat(n_row, n_col);

		SOPVector stack;
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

	SOPVector CNotGridCol(const mctdhBasis& basis, size_t n_row, size_t n_col, size_t offset, size_t shift) {
		Matrix<bool> m = true_mat(n_row, n_col);
		m.print();

		SOPVector stack;
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

	SOPVector CNotGridPattern(const mctdhBasis& basis, size_t n_row, size_t n_col, size_t which) {
		switch (which) {
			case 0:return CNotGridRow(basis, n_row, n_col, 0, 0);
			case 1:return CNotGridRow(basis, n_row, n_col, 0, 1);
			case 2:return CNotGridRow(basis, n_row, n_col, 1, 0);
			case 3:return CNotGridRow(basis, n_row, n_col, 1, 1);
			case 4:return CNotGridCol(basis, n_row, n_col, 0, 0);
			case 5:return CNotGridCol(basis, n_row, n_col, 0, 1);
			case 6:return CNotGridCol(basis, n_row, n_col, 1, 0);
			case 7:return CNotGridCol(basis, n_row, n_col, 1, 1);
			default:cerr << "Grid pattern are defined for patterns 0,..,7\n";
				exit(1);
		}
	}

	SOPVector CNotGridSwipe(const mctdhBasis& basis, size_t n_row, size_t n_col) {
		SOPVector stack;
		for (size_t n = 0; n < 8; ++n) {
			cout << "n = " << n << endl;
			SOPVector GridOp = CNotGridPattern(basis, n_row, n_col, n);
			stack.append(GridOp);
		}
		return stack;
	}

	SOPVector GoogleIterations(const mctdhBasis& basis, size_t n_row, size_t n_col, size_t n_swipes) {
		SOPVector stack;
		{
			MultiParticleOperator M;
			SOP sop;
			for (size_t n = 0; n < basis.nLeaves(); n += 2) {
				M.push_back(Hadamard, n);
			}
			sop.push_back(M, 1.0);
			stack.push_back(sop);
		}
		for (size_t n = 0; n < n_swipes; ++n) {
			auto SwipeStack = CNotGridSwipe(basis, n_row, n_col);
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
}
