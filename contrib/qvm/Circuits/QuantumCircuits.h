//
// Created by Roman Ellerbrock on 8/18/20.
//

#ifndef QUANTUMCIRCUITS_H
#define QUANTUMCIRCUITS_H

#include "../QuantumCircuit.h"
#include "../QuantumInstruction.h"
#include "../MeasurementInstruction.h"
#include "../Util/long_integer.h"
#include "../Circuits/Arithmetic.h"

namespace Circuits {
	SOPVectorcd cUa(const YAML::Node& node,
		const Register& control, const Register& x, const Register& b,
		const Register& zero, bool adjoint, size_t approx);

	SOPVectorcd cMULTmod(const YAML::Node& node,
		const Register& control, const Register& x, const Register& b,
		const Register& zero, bool adjoint);

	SOPVectorcd ccADDmod(const YAML::Node& node,
		const Register& c1, const Register& c2, const Register& b,
		const Register& zero, bool adjoint);

	shared_ptr<LeafOperatorcd> Rspecial(
		const Measurements::Measurement & m, bool adjoint);

	SOPVectorcd Rspecial(const YAML::Node& node,
		const Register& target, const QVMState& state,
		bool adjoint);

	QuantumCircuit Shor (const Register& reg,
		const cpp_int& a, const cpp_int& N);

    QuantumCircuit ShorFull (const Register& reg,
                         const cpp_int& a, const cpp_int& N);

	SOPVectorcd transversalIsing(size_t steps);
	QuantumCircuit isingModel(double p, size_t steps);
	QuantumCircuit isingIntegration(const Register& reg, size_t depth, double p, size_t steps);
}

#endif //QUANTUMCIRCUITS_H
