//
// Created by Roman Ellerbrock on 4/22/20.
//

#ifndef LABEL_H
#define LABEL_H
#include "Instruction.h"
#include "yaml-cpp/yaml.h"

class Label : public Instruction {
public:
	Label(const YAML::Node& node) {
		type_ = "Label";
		if (auto name_par = node["name"]) {
			name_ = name_par.as<string>();
		} else {
			cerr << "Set name for label.\n";
			exit(3);
		}
	}
	~Label() = default;

	string name_;

	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override {
		state.pos_++;
	}

	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override {
		state.pos_++;
	}

	void setLabels(size_t& pos, map<string, size_t>& label) override {
		label[name_] = pos;
		pos++;
	}

};


#endif //LABEL_H
