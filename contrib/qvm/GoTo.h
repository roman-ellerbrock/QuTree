//
// Created by Roman Ellerbrock on 4/23/20.
//

#ifndef GOTO_H
#define GOTO_H
#include "yaml-cpp/yaml.h"
#include "Instruction.h"

class GoTo : public Instruction {
public:
	GoTo(const YAML::Node& node) {
		type_ = "goto";
		if (auto name_par = node["name"]) {
			label_ = name_par.as<string>();
		} else {
			cerr << "Did not specify name of label for goto statement\n";
			exit(3);
		}
	}
	~GoTo() = default;

	void apply(QVMState& state, ostream& os) {
		if (!state.label_[label_]) {
			cerr << "Runtime Error: goto label does not exist: " << label_ << endl;
			exit(5);
		}
		state.pos_ = state.label_.at(label_);
	}

	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override {
		apply(state, os);
	}

	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override {
		apply(state, os);
	}


private:
	string label_;
};


#endif //GOTO_H
