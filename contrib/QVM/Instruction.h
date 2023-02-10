//
// Created by Roman Ellerbrock on 4/15/20.
//

#ifndef INSTRUCTION_H
#define INSTRUCTION_H

#include "TreeClasses/TensorTree.h"
#include "QVMState.h"
#include <map>
#include "FullRank.h"

template <typename T>
T eval(const YAML::Node& node, const string& key) {
	T val;
	if (auto par = node[key]) {
		return par.as<T>();
	} else {
		cerr << "Did not specify key '" << key << "' in yaml node " << node << endl;
		exit(3);
	}
}


/**
 * class Instruction
 *
 * \brief Instruction for a quantum computer
 *
 * This class is the instruction that is operated on a quantum computer.
 * It provides an abstract function that specifys what is to be done when
 * the instruction is applied. The class is inherited to
 * QuantumInstruction, ClassicalInstruction, GotoInstruction, GotoLabel
 *
 */

typedef map<string, string> QC_Memory;

class Instruction {
public:
	Instruction() = default;
	~Instruction() = default;

	virtual void apply(Wavefunction& Psi, MatrixTreecd& rho, QVMState& state, ostream& os) = 0;
	virtual void apply(FullRank::Wavefunction& Psi, QVMState& state, ostream& os) = 0;

	virtual void setLabels(size_t& pos, map<string, size_t>& label) { pos++; }

	string type_;
};

#endif //INSTRUCTION_H
