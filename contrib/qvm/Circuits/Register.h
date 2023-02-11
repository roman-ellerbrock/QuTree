//
// Created by Roman Ellerbrock on 2019-10-26.
//

#ifndef REGISTER_H
#define REGISTER_H
#include <string>
#include <vector>
#include <iostream>
#include <assert.h>
#include "yaml-cpp/yaml.h"
#include "Util/string_ext.h"

/**
 * \class Register
 *
 * \brief Manages register properties like sub-registers, size, start, end, etc.
 *
 */

using namespace std;
class Register {
public:
	Register()
		: size_(0), front_(0) {}

	Register(size_t begin, size_t size, string name_ = "");
	~Register() = default;

	explicit Register(const YAML::Node& node, size_t begin);

	///
	void push_back(Register reg) {
		subregs_.push_back(move(reg));
	}

	///
	size_t front() const { return front_; }

	///
	size_t size() const {
		size_t total = size_;
		for (const Register& reg : subregs_) {
			total += reg.size();
		}
		return total;
	}

	size_t largestBit()const { return front(); }
	size_t smallestBit()const { return back(); }

	void setThisSize(size_t size) { size_ = size; }

	void setBegin(size_t begin) {
		for (auto& reg: subregs_) {
			reg.setBegin(begin);
			begin += reg.size();
		}
		front_ = begin;
	}

	///
	size_t back() const { return (front_ + size() - 1); }
	size_t end() const { return (front_ + size()); }

	///
	size_t nSubRegisters() const { return subregs_.size(); }

	///
	const string& name() const { return name_; }
	string& name() { return name_; }

	/// This is a primitive Register (a leaf) if there are not subregisters
	bool isLeaf() const { return subregs_.empty(); }

/*
 Note:
 Returning a qubit as a register requires create a temporary object.
 This can either only be returned by value or stored as a member
 which would remove the const qualifier.
 The best solution I found so far was to ignore any number in the path
 and manage that afterwards from outside.
  		Register operator[](size_t i) const {
		assert(i < Size());
		auto tmp_name = name_ + "[" + to_string(i) + "]";
		return Register(i, 1, tmp_name);
	}
*/

	const Register& operator[](const string& name_)const;

	/// Can go directly to deep subregister. Ignores number inputs.
	const Register& operator[](const vector<string>& path) const {
		const Register* t = this;
		for (const auto& el : path) {
			if (is_number(el)) { break; }
			t = &t->operator[](el);
		}
		return *t;
	}

	const Register& operator[](const YAML::Node& node) const;

	Register& operator[](const string& name_);

	void print()const;

private:
	string name_;
	vector<Register> subregs_;
	shared_ptr<Register> tmp_;

	size_t front_;
	size_t size_;
};


#endif //REGISTER_H
