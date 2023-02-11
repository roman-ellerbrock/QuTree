//
// Created by Roman Ellerbrock on 4/24/20.
//

#include "IfElse.h"

ConditionArg::ConditionArg(const YAML::Node& node)
	: is_label_(true) {
	/// Could be clas1 or 0
	if (node.IsScalar()) {
		/// case label
		label_ = node.as<string>();
		if (is_number(label_)) {
			/// case number
			is_label_ = false;
			val_.push_back(node.as<size_t>());
		}
	} else {
		/// case number sequence
		is_label_ = false;
		val_ = to_vector(node);
	}
}

bool Condition::evaluate(const QVMState& state) const {
	const vector<size_t>& l = left_.eval(state);
	const vector<size_t>& r = right_.eval(state);
	assert(l.size() == r.size());
	assert(comp_ == "==" || comp_ == "!=");
	if (comp_ == "==") {
		for (size_t i = 0; i < l.size(); ++i) {
			if (l[i] != r[i]) { return false; }
		}
		return true;
	} else if (comp_ == "!=") {
		for (size_t i = 0; i < l.size(); ++i) {
			if (l[i] == r[i]) { return false; }
		}
		return true;
	} else {
		cerr << "Invalid comparator.\n";
		exit(5);
	}
}

vector<size_t> to_vector(const YAML::Node& node) {
	/// converts nodes like "[0 1 2]" or "0" into a vector<size_t>
	vector<size_t> vec;
	if (node.IsSequence()) {
		for (const auto& x : node) {
			/// Now x should be x number
			assert(x.IsScalar());
			assert(is_number(x.as<string>()));
			vec.push_back((size_t) x.as<size_t>());
		}
	} else if (node.IsScalar()) {
		assert(is_number(node.as<string>()));
		vec.push_back(node.as<size_t>());
	}
	return vec;
}

YAML::Node save_eval(const YAML::Node& node, const string& key) {
	if (node[key]) {
		return node[key];
	} else {
		cerr << "Did not find key '" << key << "' in YAML-Node:\n";
		cerr << node << endl;
		exit(2);
	}
}


ostream& operator<<(ostream& os, const ConditionArg& arg) {
	return arg.operator<<(os);
}

ostream& operator<<(ostream& os, const Condition& c) {
	return c.operator<<(os);
}

ostream& operator<<(ostream& os, const IfElse& ifelse) {
	return ifelse.operator<<(os);
}
