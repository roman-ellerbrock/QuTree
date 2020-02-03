#pragma once
#include "TreeStructuredObject.h"
#include "TensorTree.h"
#include "SingleParticleOperator.h"
#include "SparseFactorMatrixTree.h"
#include "SparseTreeStructuredObject.h"


template<typename T>
class SparseHoleMatrixTree: public SparseTreeStructuredObject<Matrix<T>>
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
{
public:
	using SparseTreeStructuredObject<Matrix<T>>::Active;
	using SparseTreeStructuredObject<Matrix<T>>::operator[];
	using SparseTreeStructuredObject<Matrix<T>>::Initialize;
	using SparseTreeStructuredObject<Matrix<T>>::attributes_;

	/// Create HoleMatrixTree from file
	SparseHoleMatrixTree(const MPO<T>& M, const TTBasis& basis, const string& filename);

	/// Create HoleMatrixTree for a given tree-marker
	SparseHoleMatrixTree(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(active_, basis) {
		Initialize(basis);
	}

	/// Create HoleMatrixTree only for relevant nodes for a given Operator
	SparseHoleMatrixTree(const MPO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(cast_to_vector_size_t(M.Modes()), basis) {
		Initialize(basis);
	}

	/// Create and calculate HoleMatrixTree
	SparseHoleMatrixTree(const TensorTree<T>& Psi, const SparseFactorMatrixTree<T>& hmat,
		const MPO<T>& M, const TTBasis& basis);

	~SparseHoleMatrixTree() = default;

	/// Create Matrices for active_ nodes in the tree
	void Initialize(const TTBasis& basis) override;

	/// Calculate Hole-Matrices from FactorMatrixTree
	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseFactorMatrixTree<T>& hmat, const TTBasis& basis);

	/// Calculate Hole-Matrices from FactorMatrixTree
	void Calculate(const TensorTree<T>& Psi,
		const SparseFactorMatrixTree<T>& hmat, const TTBasis& basis) {
		Calculate(Psi, Psi, hmat, basis);
	}

	/// Calculate Hole-Matrices from FactorMatrixTree at active nodes
	void Calculate(const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const SparseFactorMatrixTree<T>& hmat,
		const TreeMarker& act);

	/// Calculate Hole-Matrices from FactorMatrixTree at active nodes
	void Calculate(const TensorTree<T>& Psi,
		const SparseFactorMatrixTree<T>& hmat,
		const TreeMarker& act) {
		Calculate(Psi, Psi, hmat, act);
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
ostream& operator>>(ostream& os, const SparseHoleMatrixTree<T>& hmat);

template<typename T>
istream& operator<<(istream& is, SparseHoleMatrixTree<T>& hmat);

typedef SparseHoleMatrixTree<complex<double>> SparseHoleMatrixTreecd;

typedef SparseHoleMatrixTree<double> SparseHoleMatrixTreed;



