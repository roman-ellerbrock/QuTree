//
// Created by Roman Ellerbrock on 8/28/20.
//

#ifndef OUTPUTINSTRUCTION_H
#define OUTPUTINSTRUCTION_H
#include "Instruction.h"

class OutputInstruction : public Instruction {
public:
	OutputInstruction() = default;
	OutputInstruction(const YAML::Node& node, const Register& reg) {
	}
	~OutputInstruction() = default;

	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override;
	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override;

};


#endif //OUTPUTINSTRUCTION_H
