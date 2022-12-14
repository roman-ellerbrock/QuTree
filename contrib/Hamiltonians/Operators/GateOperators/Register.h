//
// Created by Roman Ellerbrock on 2019-10-26.
//

#ifndef REGISTER_H
#define REGISTER_H
#include "token.h"

/**
 * \class Register
 *
 * \brief Manages register properties like sub-registers, size, start, end, etc.
 *
 */

class Register {
public:
	Register()
		: size(0), begin(0) {}

	Register(size_t begin_, size_t size_, string name_ = "");
	~Register() = default;

	///
	void push_back(Register reg) {
		subregs.push_back(move(reg));
	}

	///
	size_t Begin() const { return begin; }

	///
	size_t Size() const {
		size_t total = size;
		for (const Register& reg : subregs) {
			total += reg.Size();
		}
		return total;
	}

	size_t LargestBit()const { return Begin(); }
	size_t SmallestBit()const { return Last(); }

	void SetThisSize(size_t size_) { size = size_; }

	void SetBegin(size_t begin_) {
		for (auto& reg: subregs) {
			reg.SetBegin(begin_);
			begin_ += reg.Size();
		}
		begin = begin_;
	}

	///
	size_t Last() const { return (begin + Size() - 1); }
	size_t End() const { return (begin + Size()); }

	///
	size_t nSubRegisters() const { return subregs.size(); }

	///
	string Name() const { return name; }

	/// This is a primitive Register (a leaf) if there are not subregisters
	bool IsLeaf() const { return subregs.empty(); }

	const Register& operator[](size_t i) const {
		assert(i < Size());
		return subregs[i];
	}

	Register& operator[](size_t i) {
		assert(i < Size());
		return subregs[i];
	}

	const Register& operator[](const string& name_)const {
		for (const Register& reg: subregs) {
			if (reg.Name() == name_) {
				return reg;
			}
		}
		cerr << "No such subregister with name: " << name_ << endl;
		exit(1);
	}

	Register& operator[](const string& name_) {
		for (Register& reg: subregs) {
			if (reg.Name() == name_) {
				return reg;
			}
		}
		cerr << "No such subregister with name: " << name_ << endl;
		exit(1);
	}

	void print()const {
		if (subregs.empty()) {
			cout << "Reg: " << name << ", begin: " << Begin() << ", last: " << Last() << ", qubits: " << Size()
				 << endl;
		} else {
			cout << "Reg: " << name << " with " << subregs.size() << " subregisters" << endl;
			for (const auto& subreg : subregs) {
				subreg.print();
			}
		}
	}

private:
	string name;
	vector<Register> subregs;

	size_t begin;
	size_t size;
};


#endif //REGISTER_H
