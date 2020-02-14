//
// Created by Roman Ellerbrock on 2/12/20.
//
#ifndef SPARSEMATRIXTREE_H
#define SPARSEMATRIXTREE_H
#include "Tree/SparseTreeStructuredObject.h"

template <typename T>
class SparseMatrixTree : public SparseTreeStructuredObject<Matrix<T>>{
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
	using SparseTreeStructuredObject<Matrix<T>>::attributes_;
public:
	using SparseTreeStructuredObject<Matrix<T>>::Active;
	using SparseTreeStructuredObject<Matrix<T>>::operator[];
	using SparseTreeStructuredObject<Matrix<T>>::Initialize;
	using SparseTreeStructuredObject<Matrix<T>>::Size;

	/// Create HoleMatrixTree for a given tree-marker
	SparseMatrixTree(shared_ptr<SubTree>& active_, const TTBasis& tree)
	: SparseTreeStructuredObject<Matrix<T>>(active_, tree) {
		Initialize(tree);
	}

	/// Create HoleMatrixTree only for relevant nodes for a given Operator
	SparseMatrixTree(const MLO<T>& M, const TTBasis& tree)
		: SparseTreeStructuredObject<Matrix<T>>(M.Modes(), tree) {
		Initialize(tree);
	}

	~SparseMatrixTree() = default;

	/// Create Matrices for active_ nodes in the tree
	void Initialize(const TTBasis& tree) override;

	/// I/O
	void print(ostream& os = cout);
	void Write(ostream& os) const;
	void Write(const string& filename) const;
	void Read(istream& is);
	void Read(const string& filename);
};

template<typename T>
ostream& operator>>(ostream& os, const SparseMatrixTree<T>& hmat);

template<typename T>
istream& operator<<(istream& is, SparseMatrixTree<T>& hmat);

typedef SparseMatrixTree<complex<double>> SparseMatrixTreecd;

typedef SparseMatrixTree<double> SparseMatrixTreed;




#endif //SPARSEMATRIXTREE_H
