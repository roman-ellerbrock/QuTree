//
// Created by Roman Ellerbrock on 5/23/20.
//

#ifndef MATRIXLISTTREE_H
#define MATRIXLISTTREE_H
#include "TreeOperators/TensorOperators/TensorOperatorTree.h"
#include "TreeClasses/TensorTree.h"

typedef vector<Matrixcd> MatrixList;

class MatrixListTree: public NodeAttribute<MatrixList> {
public:
	/*!
	 * \brief A tree graph of MatrixLists
	 *
	 * \ingroup TreeOperators
	 *
	 */

	/// Default constructor
	MatrixListTree() = default;

	/// Default Destructor
	~MatrixListTree() = default;

	/// Constructor that initializes the mlHMatrices
	MatrixListTree(const Tree& tree)
	{ initialize(tree); }

	/// Initialize the mlHMatrices
	void initialize(const Tree& tree);

	/// Print the ml-H-matrices
	void print(const Tree& tree, ostream& os = cout) const;

};


#endif //MATRIXLISTTREE_H
