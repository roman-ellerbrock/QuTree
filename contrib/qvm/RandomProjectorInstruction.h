//
// Created by Roman Ellerbrock on 1/18/21.
//

#ifndef RANDOMPROJECTORINSTRUCTION_H
#define RANDOMPROJECTORINSTRUCTION_H
#include "Instruction.h"

class RandomProjectorInstruction : public Instruction {
public:
	RandomProjectorInstruction() = default;
	~RandomProjectorInstruction() = default;

	RandomProjectorInstruction(size_t q) : q_(q) {}
	RandomProjectorInstruction(const YAML::Node& node, const Register& reg) {

	}

	void apply(Wavefunction& Psi, MatrixTreecd& rho,
		QVMState& state, ostream& os);
	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os);

private:
	size_t q_;
};


#endif //RANDOMPROJECTORINSTRUCTION_H
