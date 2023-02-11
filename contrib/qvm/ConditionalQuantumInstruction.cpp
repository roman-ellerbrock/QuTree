//
// Created by Roman Ellerbrock on 8/18/20.
//

#include "ConditionalQuantumInstruction.h"

bool is_equal(const Bitstring& x, const Bitstring& y) {
	bool res = true;
	for (size_t i = 0; i < x.size(); ++i) {
		if (x[i] != y[i]) {
			res = false;
		}
	}
	return res;
}

ConditionalQuantumInstruction::ConditionalQuantumInstruction(
	const YAML::Node& node, const Register& reg) {
	creg_name_ = eval<string>(node, "classical_reg");
	auto value = eval<string>(node, "value");
	cval_.clear();
	cout << "string from measurement:\n";
	for (char & c: value) {
		int z = atoi(&c);
		cout << z << endl;
		cval_.push_back(z);
	}
	cout << endl;

	YAML::Node n = node;
	n["instr"] = "quantum";
	qinstr_ = make_shared<QuantumInstruction>(node, reg);
}

void ConditionalQuantumInstruction::apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) {
	const Bitstring& bits = state.classical_reg_[creg_name_];

	if (is_equal(bits, cval_)) {
		qinstr_->apply(Psi, rho, state, os);
	} else {
		++state.pos_;
	}
}

void ConditionalQuantumInstruction::apply(FullRank::Wavefunction& Psi,
	QVMState& state, ostream& os) {
	const Bitstring& bits = state.classical_reg_[creg_name_];
	if (is_equal(bits, cval_)) {
		qinstr_->apply(Psi, state, os);
	} else {
		++state.pos_;
	}
}

ConditionalQuantumInstruction::ConditionalQuantumInstruction(const string& classical_reg,
	const string& value, const string& circuit_name,
	const vector<Register>& targets) : creg_name_(classical_reg){

	for (char c: value) {
		int z = atoi(&c);
		cout << z << endl;
		cval_.push_back(z);
	}
	qinstr_ = make_shared<QuantumInstruction>(circuit_name, targets);
}
