//
// Created by Roman Ellerbrock on 2/27/20.
//

#ifndef MCTDH_STATE_H
#define MCTDH_STATE_H
#include "TreeShape/Tree.h"
#include "TreeOperators/Hamiltonian.h"

typedef TensorTreecd Wavefunction;

struct mctdh_state {
	mctdh_state():rng_(0) {}
//	mctdh_state():rng_(time(nullptr)) {}
	~mctdh_state() = default;

	Tree tree_;
	Tree cdvrtree_;
	shared_ptr<Hamiltonian> hamiltonian_;
	map<string, Wavefunction> wavefunctions_;

	mt19937 rng_;

	void print(ostream& os = cout) const {
		os << "Tree: " << endl;
		tree_.print(os);
		os << "Hamiltonian:" << endl;
		hamiltonian_->print(os);
		os << "Wavefunctions:" << endl;
		for (const auto& pair : wavefunctions_) {
			os << pair.first << endl;
			pair.second.print(tree_, os);
		}
	}

};


#endif //MCTDH_STATE_H
