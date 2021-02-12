//
// Created by Roman Ellerbrock on 2019-07-13.
//

#include "TreeOperators/TreeStructured/TreeSOP.h"

void NodeOperator::Print(const vector<string>& names, size_t indent) const {
	for (size_t i = 0; i < indent; ++i) {
		cout << "\t";
	}
	if (names.empty()) {
		cout << "(" << part() << ", " << mode() << ")";
	} else {
		cout << "(" << names[part()] << ", " << mode() << ")";
	}
}

void print(const NodeProductOperator& mpo, complex<double> coeff, const vector<string>& names, size_t indent) {
	// print coeff
	if (imag(coeff) == 0.) {
		cout << "\t" << real(coeff) << " * ";
	} else {
		cout << "\t" << coeff << " * ";
	}
	// print spos
	for (size_t i = 0; i < mpo.size(); ++i) {
		if (i > 0) { cout << " * "; }
		const NodeOperator& spo = mpo[i];
		spo.Print(names, indent);
	}
//		cout << "\n";
}

void print(const NodeSOP& sop, const vector<string>& names, size_t indent) {
	for (size_t i = 0; i < indent; ++i) {
		cout << "\t";
	}
	for (size_t l = 0; l < sop.size(); ++l) {
//		cout << "l=" << l << ":";
		if (l > 0) { cout << " + "; } else { cout << "   "; }
		::print(sop[l], sop.Coeff(l), names, indent);
	}
	cout << "\n";
}

void print(const NodeSOPlist& sopl, const vector<string>& names, size_t indent) {
	for (size_t n = 0; n < sopl.size(); ++n) {
		cout << "SOP No. " << n << ":\n";
		const NodeSOP& sop = sopl[n];
		::print(sop, names, 0);
		cout << "\n";
	}
	cout << "\n";
}

void TreeSOP::print(const Tree& tree) const {
	cout << "Printing multilayer SOP operator:\n";
	for (const Node& node : tree) {
		node.info();
		if (node.isBottomlayer()) {
			const NodeSOPlist& sopl = operator[](node);
			const Leaf& leaf = node.getLeaf();
			::print(sopl, leafoperatornames[leaf.Mode()], 0);
		} else {
			const NodeSOPlist& sopl = operator[](node);
			::print(sopl, vector<string>(), 0);
		}
	}

	cout << "Bottomlayer operator library:\n";
	size_t mode = 0;
	for (const auto& names : leafoperatornames) {
		cout << "mode = " << mode << "\t";
		for (const auto& name : names) {
			cout << name << "\t";
		}
		cout << endl;
		mode++;
	}
}

bool IsActive(const NodeOperator& h, const NodeProductOperator& M) {
	for (const NodeOperator& hm : M) {
		if (h == hm) { return true; }
	}
	return false;
}


