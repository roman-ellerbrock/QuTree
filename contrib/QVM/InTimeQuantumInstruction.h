//
// Created by Roman Ellerbrock on 8/20/20.
//

#ifndef INTIMEQUANTUMINSTRUCTION_H
#define INTIMEQUANTUMINSTRUCTION_H
#include "QuantumInstruction.h"

class InTimeQuantumInstruction : public Instruction {
public:
	InTimeQuantumInstruction(const YAML::Node& node, const Register& reg);
	~InTimeQuantumInstruction() = default;

	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override;
	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override;

private:
	QuantumInstruction qinstr_;
	YAML::Node node_;
};


#endif //INTIMEQUANTUMINSTRUCTION_H
