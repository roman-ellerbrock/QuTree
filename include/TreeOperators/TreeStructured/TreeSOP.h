//
// Created by Roman Ellerbrock on 5/21/20.
//

#ifndef TREESOP_H
#define TREESOP_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeOperators/LeafOperator.h"

/**
 * \class TreeSOP
 * \ingroup TreeOperators
 * \brief A Tree-structured Sum-of-Product (SOP) Operator.
 *
 * */

typedef vector<shared_ptr<LeafOperatorcd>> LeafOperatorLib;

//typedef pair<size_t, size_t> NodeOperator;

class NodeOperator {
public:
	NodeOperator(size_t part, size_t mode)
		: mode_(mode), part_(part) {}

	~NodeOperator() = default;

	size_t mode() const { return mode_; }

	size_t part() const { return part_; }

	bool operator==(const NodeOperator& h) const {
		return ((h.part() == part()) && (h.mode() == mode()));
	}

	void Print(const vector<string>& names, size_t indent = 0) const;

private:
	size_t mode_;
	size_t part_;
};

typedef vector<NodeOperator> NodeProductOperator;

bool IsActive(const NodeOperator& h, const NodeProductOperator& M);

class NodeSOP: public vector<NodeProductOperator> {
public:
	NodeSOP() = default;
	~NodeSOP() = default;

	void push_back(const NodeProductOperator& M, const complex<double> coeff) {
		vector<NodeProductOperator>::push_back(M);
		coeffs.push_back(coeff);
	}

	complex<double> Coeff(size_t part) const {
		assert(part < coeffs.size());
		return coeffs[part];
	}

private:
	vector<complex<double>> coeffs;
};

typedef vector<NodeSOP> NodeSOPlist;

// Output routines
void print(const NodeSOPlist& sopl, const vector<string>& names, size_t indent);
void print(const NodeSOP& sop, const vector<string>& names, size_t indent);

class TreeSOP
	: public NodeAttribute<NodeSOPlist> {
public:
	TreeSOP() = default;

	explicit TreeSOP(const Tree& tree) {
		Initialize(tree);
	}

	~TreeSOP() = default;

	void Initialize(const Tree& tree) {
		// Clear Bottomlayer Operators and initialize them
		leafoperatorlibs.clear();
		leafoperatornames.clear();
		SpecialInitializeBottom(tree);

		// Clear Upper layer_ Operators and initialize them
		attributes_.resize(tree.nNodes());
		SpecialInitialize(tree);
	}

	LeafOperatorLib& leafOperatorLib(const Node& node) {
		assert(node.isBottomlayer());
		const Leaf& phys = node.getLeaf();
		const size_t mode = phys.mode();
		assert(mode < leafoperatorlibs.size());
		return leafoperatorlibs[mode];
	}

	const LeafOperatorLib& leafOperatorLib(const Node& node) const {
		assert(node.isBottomlayer());
		const Leaf& leaf = node.getLeaf();
		const size_t mode = leaf.mode();
		assert(mode < leafoperatorlibs.size());
		return leafoperatorlibs[mode];
	}

	void push_back(const LeafOperatorLib& lib,
		const vector<string>& names) {
		assert(lib.size() == names.size());
		leafoperatorlibs.push_back(lib);
		leafoperatornames.push_back(names);
	}

	void print(const Tree& tree) const;

	vector<vector<string>> LeafOperatorLibNames() const {
			return leafoperatornames;
	}

	bool empty()const { return attributes_.empty(); }

private:
	vector<LeafOperatorLib> leafoperatorlibs;

	vector<vector<string>> leafoperatornames;

	virtual void SpecialInitializeBottom(const Tree& tree) {};

	virtual void SpecialInitialize(const Tree& tree) {};
};

#endif //TREESOP_H
