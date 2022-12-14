//
// Created by Roman Ellerbrock on 2019-05-10.
//
#include "QFT.h"

namespace GateQFT {
	SOPVector QFT(const Register& reg, bool adjungate, size_t approx) {
		return QFT(reg.Begin(), reg.Size(), adjungate, approx);
	}

	SOPVector iQFT(const Register& reg, bool adjungate, size_t approx) {
		return iQFT(reg.Begin(), reg.Size(), adjungate, approx);
	}

	sopList ControlledRotation(const Register& reg1, const Register& reg2,
		bool adjungate, size_t approx) {
		assert(reg1.Size() == reg2.Size());
		return ControlledRotation(reg1.Size(), reg1.Begin(), reg2.Begin(), adjungate, approx);
	}

	sopList QFT(size_t mode_start, size_t n_bit, bool adjungate, size_t approx) {
		using namespace GateOperators;
		sopList qft;
		size_t mode_end = mode_start + n_bit;
		for (size_t n = 0; n < n_bit; ++n) {
			size_t bit_a = mode_start  + n;
			{
				MPO M(Hadamard, bit_a);
				SOP sop(M, 1.);
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
	}

	sopList iQFT(size_t mode_start, size_t n_bit, bool adjungate, size_t approx) {
		auto iqft = QFT(mode_start, n_bit, !adjungate, approx);
		std::reverse(iqft.begin(), iqft.end());
		return iqft;
	}

	sopList ControlledRotation(size_t n_bit, size_t mode_a, size_t mode_b, bool adjungate, size_t approx) {
		sopList stack;

		for (size_t b = 0; b < n_bit; ++b) {
			for (size_t a = b; a < n_bit; ++a) {
				size_t na = mode_a + a;
				size_t nb = mode_b + b;
				size_t angle_k = a - b + 1;
				if (approx == 0 || angle_k < approx) {
					stack.emplace_back(GateOperators::Uk(na, nb, angle_k, adjungate));
//					cout << "U(" << na << ", " << nb << ")[" << angle_k << "] " << endl;
				}
			}
		}
		return stack;
	}
}


