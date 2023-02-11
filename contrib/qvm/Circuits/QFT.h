//
// Created by Roman Ellerbrock on 2019-05-10.
//

#ifndef MCTDH_QFT_H
#define MCTDH_QFT_H
#include "GateOperators.h"
#include "TreeOperators/SOPVector.h"

namespace Circuits {
	SOPVectorcd QFT(const Register& reg, bool adjungate = false, size_t approx = 0);
	SOPVectorcd iQFT(const Register& reg, bool adjungate = false, size_t approx = 0);
	SOPVectorcd controlledRotation(const Register& a, const Register& b, bool adjungate, size_t approx = 0);

	SOPVectorcd QFT(size_t mode_start, size_t n_bit, bool adjungate = false, size_t approx = 0);
	SOPVectorcd iQFT(size_t mode_start, size_t n_bit, bool adjungate = false, size_t approx = 0);
	SOPVectorcd controlledRotation(size_t n_bit, size_t mode_a, size_t mode_b, bool adjungate, size_t approx = 0);
}

#endif //MCTDH_QFT_H
