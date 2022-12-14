//
// Created by Roman Ellerbrock on 2019-05-10.
//

#ifndef MCTDH_QFT_H
#define MCTDH_QFT_H

#include "GateOperators.h"

namespace GateQFT {
	sopList QFT(const Register& reg, bool adjungate = false, size_t approx = 0);
	sopList iQFT(const Register& reg, bool adjungate = false, size_t approx = 0);
	sopList ControlledRotation(const Register& a, const Register& b, bool adjungate, size_t approx = 0);

	sopList QFT(size_t mode_start, size_t n_bit, bool adjungate = false, size_t approx = 0);
	sopList iQFT(size_t mode_start, size_t n_bit, bool adjungate = false, size_t approx = 0);
	sopList ControlledRotation(size_t n_bit, size_t mode_a, size_t mode_b, bool adjungate, size_t approx = 0);
}

#endif //MCTDH_QFT_H
