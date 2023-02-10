//
// Created by Roman Ellerbrock on 9/4/20.
//

#ifndef ENTROPYTREE_H
#define ENTROPYTREE_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeShape/Tree.h"
#include "TreeClasses/TensorTree.h"

class EntropyTree : public NodeAttribute<double> {
public:
	EntropyTree() = default;
	~EntropyTree() = default;
	EntropyTree(const Tree& tree) {
		initialize(tree);
	}

	void initialize(const Tree& tree);
	void calculate(const TensorTreecd& Psi, const Tree& tree);
	void perplexity(const TensorTreecd& Psi, const Tree& tree,
		size_t power = 1);

	void print(const Tree& tree);

	double metric(const Tree& tree) const;
};


#endif //ENTROPYTREE_H
