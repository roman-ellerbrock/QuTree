//
// Created by Roman Ellerbrock on 3/18/20.
//

#ifndef DELTAVTREE_H
#define DELTAVTREE_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeShape/Tree.h"

class DeltaVTree : public NodeAttribute<Tensorcd>{
public:
	explicit DeltaVTree(const Tree& tree) { Initialize(tree); }
	~DeltaVTree() = default;

	void Initialize(const Tree& tree);

	void print(const Tree& tree);

};


#endif //DELTAVTREE_H
