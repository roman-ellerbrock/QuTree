//
// Created by Roman Ellerbrock on 8/20/20.
//

#include "InTimeQuantumInstruction.h"


InTimeQuantumInstruction::InTimeQuantumInstruction(
	const YAML::Node& node, const Register& reg)
	: node_(node) { }

void InTimeQuantumInstruction::apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) {
	qinstr_ = QuantumInstruction(node_, state);
	qinstr_.apply(Psi, rho, state, os);
}

void InTimeQuantumInstruction::apply(FullRank::Wavefunction& Psi,
	QVMState& state, ostream& os) {
	qinstr_ = QuantumInstruction(node_, state);
	qinstr_.apply(Psi, state, os);
}


