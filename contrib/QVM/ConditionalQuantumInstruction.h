//
// Created by Roman Ellerbrock on 8/18/20.
//

#ifndef CONDITIONALQUANTUMINSTRUCTION_H
#define CONDITIONALQUANTUMINSTRUCTION_H
#include "QuantumInstruction.h"

class ConditionalQuantumInstruction : public Instruction {
public:
	ConditionalQuantumInstruction(const YAML::Node& node, const Register& reg);
	ConditionalQuantumInstruction(const string& classical_reg, const string& value,
		const string& circuit_name, const vector<Register>& targets);
	~ConditionalQuantumInstruction() = default;

	///
	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override;

	///
	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override;

protected:
	shared_ptr<QuantumInstruction> qinstr_;
//	QuantumInstruction qinstr_;
	string creg_name_;
	Bitstring cval_;
};


#endif //CONDITIONALQUANTUMINSTRUCTION_H
