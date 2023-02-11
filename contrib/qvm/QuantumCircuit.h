//
// Created by Roman Ellerbrock on 4/21/20.
//

#ifndef QUANTUMCIRCUIT_H
#define QUANTUMCIRCUIT_H
#include "yaml-cpp/yaml.h"
#include "Instruction.h"
#include "Circuits/Register.h"

class QuantumCircuit : public vector<shared_ptr<Instruction>> {
public:
	QuantumCircuit() = default;
	~QuantumCircuit() = default;

	QuantumCircuit(const YAML::Node& node, const Register& reg);

	virtual void execute(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state) const;

	virtual void execute(FullRank::Wavefunction& Psi, QVMState& state) const;

	string name_;
};

void run(FullRank::Wavefunction& Psi, QVMState& state, const QuantumCircuit& prog);
void run(Wavefunction& Psi, QVMState& state, const QuantumCircuit& prog);


#endif //QUANTUMCIRCUIT_H
