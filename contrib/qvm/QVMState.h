//
// Created by Roman Ellerbrock on 4/22/20.
//

#ifndef QVMSTATE_H
#define QVMSTATE_H
#include "Circuits/Register.h"
#include "yaml-cpp/yaml.h"
#include "QVMParameters.h"
#include "Applications/TreeApplyOperatorDynamic.h"
#include "FullRank.h"

typedef vector<size_t> Bitstring;

struct QVMState {
	QVMState() : pos_(0), gen_(time(NULL)) {
	}

	size_t pos_;

	Register reg_;
	map<string, Bitstring> classical_reg_;
	map<string, Wavefunction> wavefunction_;
	map<string, MatrixTreecd> rho_;
	map<string, size_t> label_;
	Tree tree_;
	Tree fr_tree_;
	Fidelity f_;

	FullRank::Wavefunction frpsi_;

	QVMParameters par_;

	mt19937 gen_;
};


#endif //QVMSTATE_H
