//
// Created by Grace Johnson on 4/21/20.
//

#ifndef QUTREEVM_RANDOM_H
#define QUTREEVM_RANDOM_H
#include "GateOperators.h"
#include "TreeOperators/SOPVector.h"
#include "../QuantumCircuit.h"
#include "../QuantumInstruction.h"
#include "../MeasurementInstruction.h"

namespace Circuits {
    SOPVectorcd random1D(const Register& reg, size_t depth, bool adjungate);

    QuantumCircuit random1D(const Register& reg, size_t depth, double p, size_t seed);

	shared_ptr<MeasurementInstruction> randomMeasurement(mt19937& gen,
		const Register& reg, double p);

	SOPVectorcd random2D(const Register& reg, size_t depth, size_t seed, bool adjoint);
	QuantumCircuit random2D(const Register& reg, size_t depth, mt19937& gen, double p);

	QuantumCircuit random2D(const Register& reg, size_t depth, mt19937& gen, double p);

	/// generate 2^num random numbers
	SOPVectorcd quantumRNG(const Register& reg, size_t num, size_t depth, mt19937& gen);

//	template <typename ... Args>
//	shared_ptr<QuantumInstruction> translate(const function<SOPVectorcd(Args ..., bool)>& f, Args ... args);
}

#endif //QUTREEVM_RANDOM_H
