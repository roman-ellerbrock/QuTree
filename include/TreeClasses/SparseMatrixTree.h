//
// Created by Roman Ellerbrock on 2/12/20.
//
#ifndef SPARSEMATRIXTREE_H
#define SPARSEMATRIXTREE_H
#include "TreeClasses/SparseNodeAttribute.h"
#include "TreeClasses/TensorTree.h"

template <typename T>
class SparseMatrixTree : public SparseNodeAttribute<Matrix<T>>{
/**
 * \class MatrixTree
 *
 * \ingroup Tree-Classes
 *
 * \brief This class represents a tree with matrices asigned to every node.
 *
 * The Hole-Matrices are the result of hole-products of tensor trees.
 * In a physical context, the hole-matrices are representation of
 * mean-field operators when working with tensor tree wavefunctions.
 * */
	using SparseNodeAttribute<Matrix<T>>::attributes_;
public:
	using SparseNodeAttribute<Matrix<T>>::sparseTree;
	using SparseNodeAttribute<Matrix<T>>::operator[];
	using SparseNodeAttribute<Matrix<T>>::initialize;
	using SparseNodeAttribute<Matrix<T>>::size;

	/// Create HoleMatrixTree for a given tree-marker
	SparseMatrixTree(shared_ptr<SparseTree> active_, const Tree& tree)
	: SparseNodeAttribute<Matrix<T>>(active_, tree) {
		initialize(tree);
	}

	explicit SparseMatrixTree(const Tree& tree)
		: SparseMatrixTree(make_shared<SparseTree>(tree), tree) {}

	/// Create HoleMatrixTree only for relevant nodes for a given Operator
	SparseMatrixTree(const MLO<T>& M, const Tree& tree, bool tail = true, bool inverse_tree = false)
		: SparseNodeAttribute<Matrix<T>>(M.targetLeaves(), tree, tail, inverse_tree) {
		initialize(tree);
	}

	/// Create HoleMatrixTree only for relevant nodes for a given Operator
	SparseMatrixTree(const vector<size_t>& targets, const Tree& tree, bool tail = true, bool inverse_tree = false)
		: SparseNodeAttribute<Matrix<T>>(targets, tree, tail, inverse_tree) {
		initialize(tree);
	}

	~SparseMatrixTree() = default;

	/// Create Matrices for active_ nodes in the tree
	void initialize(const Tree& tree) override;

	/// I/O
	void print(ostream& os = cout) const;
	void write(ostream& os) const;
	void write(const string& filename) const;
	void read(istream& is);
	void read(const string& filename);
};

template<typename T>
ostream& operator>>(ostream& os, const SparseMatrixTree<T>& hmat);

template<typename T>
istream& operator<<(istream& is, SparseMatrixTree<T>& hmat);

typedef SparseMatrixTree<complex<double>> SparseMatrixTreecd;
typedef SparseMatrixTree<double> SparseMatrixTreed;

template <typename T>
using SparseMatrixTrees = vector<SparseMatrixTree<T>>;

typedef SparseMatrixTrees<complex<double>> SparseMatrixTreescd;
typedef SparseMatrixTrees<double> SparseMatrixTreesd;



#endif //SPARSEMATRIXTREE_H
