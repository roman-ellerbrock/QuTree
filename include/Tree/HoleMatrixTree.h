#pragma once
#include "TreeStructuredObject.h"
#include "TensorTree.h"
#include "SingleParticleOperator.h"
#include "FactorMatrixTree.h"
#include "SparseTreeStructuredObject.h"

/**
 * \class HoleMatrixTree
 *
 * \ingroup Tree-Classes
 *
 * \brief This class represents the Hole-Matrices for Trees.
 *
 * The Hole-Matrices are the result of hole-products of tensor trees.
 * In a physical context, the hole-matrices are representation of
 * mean-field operators when working with tensor tree wavefunctions.
 * */

template<typename T>
class HoleMatrixTree: public SparseTreeStructuredObject<Matrix<T>> {
public:
	using SparseTreeStructuredObject<Matrix<T>>::Active;
	using SparseTreeStructuredObject<Matrix<T>>::operator[];
	using SparseTreeStructuredObject<Matrix<T>>::Initialize;
	using SparseTreeStructuredObject<Matrix<T>>::attributes;

	/// Create HoleMatrixTree from file
	HoleMatrixTree(const MPO<T>& M, const TTBasis& basis, const string& filename);

	/// Create HoleMatrixTree for a given tree-marker
	HoleMatrixTree(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(active_, basis) {
		Initialize(basis);
	}

	/// Create HoleMatrixTree only for relevant nodes for a given Operator
	HoleMatrixTree(const MPO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(cast_to_vector_size_t(M.Modes()), basis) {
		Initialize(basis);
	}

	/// Create and calculate HoleMatrixTree
	HoleMatrixTree(const TensorTree<T>& Psi, const FactorMatrixTree<T>& hmat,
		const MPO<T>& M, const TTBasis& basis);

	~HoleMatrixTree() = default;

	/// Create Matrices for active nodes in the tree
	void Initialize(const TTBasis& basis) override;

	/// Calculate Hole-Matrices form FactorMatrixTree
	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const FactorMatrixTree<T>& hmat, const TTBasis& basis);

	/// Calculate Hole-Matrices form FactorMatrixTree
	void Calculate(const TensorTree<T>& Psi,
		const FactorMatrixTree<T>& hmat, const TTBasis& basis) {
		Calculate(Psi, Psi, hmat, basis);
	}

	/// Apply HoleMatrix locally to a Tensor
	Tensor<T> Apply(const Tensor<T>& Phi, const Node& node) const;

	/// I/O
	void print(ostream& os = cout);
	void Write(ostream& os) const;
	void Write(const string& filename) const;
	void Read(istream& is);
	void Read(const string& filename);
};

template<typename T>
ostream& operator>>(ostream& os, const HoleMatrixTree<T>& hmat);

template<typename T>
istream& operator<<(istream& is, HoleMatrixTree<T>& hmat);

typedef HoleMatrixTree<complex<double>> HoleMatrixTreecd;

typedef HoleMatrixTree<double> HoleMatrixTreed;



