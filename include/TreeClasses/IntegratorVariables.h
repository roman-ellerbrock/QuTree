//
// Created by Roman Ellerbrock on 3/8/20.
//

#ifndef INTEGRATORVARIABLES_H
#define INTEGRATORVARIABLES_H
#include "TreeOperators/Hamiltonian.h"

struct IntegratorVariables {
	IntegratorVariables(double time_now_, double time_end_, double dt_, double out_,
		double accuracy_root_, double accuracy_leaf_, Wavefunction& Psi_,
		const Hamiltonian& h_, const Tree& tree_, const Tree& cdvrtree_,
		string ofname_, string ifname_, bool append,
		bool cmf_bottom = true, bool cmf_upper = true, bool cmf_top = true)
		: time_now(time_now_), time_end(time_end_), dt(dt_), out(out_),
		accuracy_root(accuracy_root_), accuracy_leaf(accuracy_leaf_), psi(&Psi_),
		h(&h_), tree(&tree_), cdvrtree(&cdvrtree_), ofname(move(ofname_)),
		ifname(move(ifname_)), append_(append),
		cmf_bottom_(cmf_bottom), cmf_upper_(cmf_upper), cmf_top_(cmf_top)
		{}

	double time_now;
	double time_end;
	double dt;
	double out;

	double accuracy_root;
	double accuracy_leaf;

	TensorTreecd* psi;
	const Hamiltonian* h;
	const Tree* tree;
	const Tree* cdvrtree;

	string ofname;
	string ifname;

	bool append_;

	bool cmf_bottom_;
	bool cmf_upper_;
	bool cmf_top_;
};


#endif //INTEGRATORVARIABLES_H
