//
// Created by Roman Ellerbrock on 4/24/20.
//

#ifndef IFELSE_H
#define IFELSE_H
#include "Instruction.h"
#include "yaml-cpp/yaml.h"
#include "QVMState.h"

vector<size_t> to_vector(const YAML::Node& node);

YAML::Node save_eval(const YAML::Node& node, const string& key);

class ConditionArg {
public:
	ConditionArg() : is_label_(true) {}

	/// \brief left- or right- argument in a condition. Can be a value or a classical reg. label
	explicit ConditionArg(const YAML::Node& node);

	const vector<size_t>& eval(const QVMState& state) const {
		if (is_label_) {
			return state.classical_reg_.at(label_);
		} else {
			return val_;
		}
	}

	ostream& operator<<(ostream& os) const {
		if (is_label_) {
			os << label_;
		} else {
			os << "[";
			for (const auto& x : val_) {
				os << " " << x;
			}
			os << "]";
		}
		return os;
	}

private:
	vector<size_t> val_;
	string label_;
	bool is_label_;
};

ostream& operator<<(ostream& os, const ConditionArg& arg);

class Condition {
public:
	explicit Condition(const YAML::Node& node) {
		assert(node["left"]);
		left_ = ConditionArg(node["left"]);

		assert(node["right"]);
		right_ = ConditionArg(node["right"]);

		comp_ = eval<string>(node, "comp");
		assert(comp_ == "==" || comp_ == "!=");
	}

	bool evaluate(const QVMState& state) const;

	ostream& operator<<(ostream& os) const {
		os << "(" << left_ << " " << comp_ << " " << right_ << ")";
		return os;
	}

	ConditionArg left_;
	ConditionArg right_;
	string comp_;
};

ostream& operator<<(ostream& os, const Condition& c);

class IfElse: public Instruction {
public:
	IfElse(const YAML::Node& node)
		: condition_(save_eval(node, "condition")) {
		type_ = "ifelse";
		if_label_ = eval<string>(node, "if");
		else_label_ = eval<string>(node, "else");
	}

	~IfElse() = default;

	void Apply(QVMState& state, ostream& os) {
		string label;

		operator<<(cout);

		if (condition_.evaluate(state)) {
			label = if_label_;
		} else {
			label = else_label_;
		}

		if (!state.label_[label]) {
			cerr << "Runtime Error: goto label does not exist: " << label << endl;
			exit(5);
		}
		state.pos_ = state.label_.at(label);
	}

	void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) override {
		Apply(state, os);
	}

	void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) override {
		Apply(state, os);
	}

	void setLabels(size_t& pos, map<string, size_t>& label) override {
		pos++;
	}

	ostream& operator<<(ostream& os) const {
		os << "if " << condition_ << " {" << endl;
		os << "\t goto: " << if_label_ << endl;
		os << "} else {" << endl;
		os << "\t goto: " << else_label_ << endl;
		os << "}" << endl;
		return os;
	}

private:
	Condition condition_;
	string if_label_, else_label_;
};

ostream& operator<<(ostream& os, const IfElse& ifelse);

#endif //IFELSE_H
