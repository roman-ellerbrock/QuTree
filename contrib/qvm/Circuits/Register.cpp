//
// Created by Roman Ellerbrock on 2019-10-26.
//

#include "Register.h"

Register::Register(size_t begin__, size_t size__, string name)
: front_(begin__), size_(size__), name_(move(name)){
}

Register::Register(const YAML::Node& node, size_t begin) {
	front_ = 0;
	string name;
	if (auto name_par = node["name"]) {
		name = name_par.as<string>();
	}

	if (node["size"]) {
		/// Primitive register
		if (node["begin"]) {
			begin = node["begin"].as<size_t>();
		}
		assert(node["size"]);
		auto size = node["size"].as<size_t>();
		*this = Register(begin, size, name);
	} else if (node["sub"]){
		/// Non-primitive reg
		assert(node.IsMap());
		Register reg;
		reg.name_ = name;
		size_t next = begin;
		for (const YAML::Node& r : node["sub"]) {
			reg.push_back(Register(r, next));
			next = reg.back() + 1;
		}
		*this = reg;
	} else {
		cerr << "Neither primitive register nor super-register.\n";
		exit(12);
	}
}

Register& Register::operator[](const string& name) {
	/// Don't ask for a number here.
	assert(!is_number(name));

	/// First check this reg, then subregs
	if (name_ == name) {
		return *this;
	}
	for (Register& reg: subregs_) {
		if (reg.name() == name) {
			return reg;
		}
	}
	cerr << "No such subregister with name: " << name << endl;
	exit(1);
}

const Register& Register::operator[](const string& name)const {
	/// If asking for a number, return the n-th element in the reg
	if (is_number(name)) {
		size_t idx = stoi(name);
		return (*this);
//		return this->operator[](idx);
	}

	/// First check this reg, then subregs
	if (name_ == name) {
		return *this;
	}
	for (const Register& reg: subregs_) {
		if (reg.name() == name) {
			return reg;
		}
	}
	cerr << "No such (sub)register with name: " << name << endl;
	exit(1);
}

const Register& Register::operator[](const YAML::Node& node) const {
	if (node.IsScalar()) {
		auto reg_name = node.as<string>();
		if (reg_name == name_) {
			return *this;
		}
		for (const auto &subreg : subregs_) {
			if (subreg.name_ == reg_name) {
				return subreg;
			}
		}
		cerr << "Cannot find register\n";
		exit(3);
	}
	if (node.IsSequence()) {
		vector<string> path;
		for (auto x : node) {
			if (!is_number(x.as<string>())) {
				path.push_back(x.as<string>());
			}
		}
		return this->operator[](path);
	}
	cerr << "Cannot process this node type\n";
	exit(3);
}

void Register::print()const {
	if (subregs_.empty()) {
		cout << "Reg: " << name_ << ", front: " << front() << ", back: " << back() << ", qubits: " << size()
			 << endl;
	} else {
		cout << "Reg: " << name_ << " with " << subregs_.size() << " subregisters" << endl;
		for (const auto& subreg : subregs_) {
			subreg.print();
		}
	}
}

