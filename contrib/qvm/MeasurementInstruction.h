//
// Created by Roman Ellerbrock on 4/23/20.
//

#ifndef MEASUREMENTINSTRUCTION_H
#define MEASUREMENTINSTRUCTION_H
#include "Instruction.h"
#include "yaml-cpp/yaml.h"

class MeasurementInstruction: public Instruction {
public:
	MeasurementInstruction(const YAML::Node& node, const Register& reg);
	MeasurementInstruction(const Register& qreg, const string& creg, bool append = false);
	MeasurementInstruction(const vector<size_t>& targets, const string& creg, bool append = false);

	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override;
	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override;

private:
	bool append_;
	string classical_reg_name_;
	vector<size_t> targets_;
};


#endif //MEASUREMENTINSTRUCTION_H
