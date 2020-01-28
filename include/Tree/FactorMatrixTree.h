/**
 * \class HMatrices
 *
 * \ingroup MCTDH-Matrices
 *
 * \brief This class represents the H-Matrices.
 *
 * The H-Matrices are objects that occur in the
 * equations of motion (EOM) of the MCTDH approach. A H-Matrix
 * is the representation of a MultiParticleOperator in the SPF-basis
 * of a given TensorTree<T>.
 * */
#pragma once
#include "Core/Matrix.h"
#include "MultiParticleOperator.h"
#include "TensorTreeBasis/TensorTreeBasis.h"
#include "TensorTree.h"
#include "SparseTreeStructuredObject.h"
#include "Core/FactorMatrix.h"

/**
 * \class FactorMatrixTree
 *
 * \ingroup Tree-Classes
 *
 * \brief This class represents the FactorMatrices for Trees.
 *
 * The FactorMatrices are the result of Tensor-products of tensor trees.
 * In a physical context, the hole-matrices are representation of
 * mean-field operators when working with tensor tree wavefunctions.
 * */

vector<size_t> cast_to_vector_size_t(const vector<int>& a);

template<typename T>
class FactorMatrixTree: public SparseTreeStructuredObject<FactorMatrix<T>> {
public:
	using SparseTreeStructuredObject<FactorMatrix<T>>::Active;
	using SparseTreeStructuredObject<FactorMatrix<T>>::operator[];
	using SparseTreeStructuredObject<FactorMatrix<T>>::Initialize;
	using SparseTreeStructuredObject<FactorMatrix<T>>::attributes;

	/// Create FactorMatrices for relevant nodes when representing an operator
	FactorMatrixTree(const MPO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<FactorMatrix<T>>(cast_to_vector_size_t(M.Modes()), basis) {
		Initialize(basis);
	}

	/// Create and calculate FactorMatrixTree for an operator
	FactorMatrixTree(const TensorTree<T>& Psi, const MPO<T>& M, const TTBasis& basis)
		: FactorMatrixTree(M, basis) {
		Calculate(Psi, M, basis);
	}

	/// Create FactorMatrices for externally marked nodes
	FactorMatrixTree(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<FactorMatrix<T>>(active_, basis) {
		Initialize(basis);
	}

	/// Read FactorMatrixTree from file
	FactorMatrixTree(const MPO<T>& M, const TTBasis& basis, const string& filename)
		: FactorMatrixTree(M, basis) {
		Read(filename);
	}

	~FactorMatrixTree() = default;

	/// Fill FactorMatrices at internally marked nodes. TreeMarker must be initialized
	void Initialize(const TTBasis& basis) override;

	/// Calculate (cross-)TreeMatrix-representation of an operator
	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const MPO<T>& M, const TTBasis& basis);

	/// Calculate TreeMatrix-representation of an operator
	void Calculate(const TensorTree<T>& Psi,
		const MPO<T>& M, const TTBasis& basis) {
		Calculate(Psi, Psi, M, basis);
	}

	/// Calculate matrices locally at a node
	void CalculateLayer(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const MPO<T>& M, const Node& node);

	/// Apply factor matrices locally
	Tensor<T> Apply(const Tensor<T>& Phi, const MPO<T>& M, const Node& node) const;

	Tensor<T> ApplyUpper(Tensor<T> Phi, const Node& node) const;

	Tensor<T> ApplyHole(Tensor<T> Phi, const Node& hole_node) const;

	/// I/O functions
	/// print human readable
	void print(const TTBasis& basis, ostream& os = cout) const;
	/// Write in readable format
	void Write(ostream& os) const;
	void Write(const string& filename) const;
	/// Read previously written (Write(..)) FactorMatrixTree
	void Read(istream& is);
	void Read(const string& filename);

protected:
	FactorMatrix<T> CalculateUpper(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const Node& node);

	FactorMatrix<T> CalculateBottom(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const MPO<T>& M, const Node& node, const Leaf& phys);
};

template<typename T>
ostream& operator>>(ostream& os, const FactorMatrixTree<T>& hmat);

template<typename T>
istream& operator<<(istream& is, FactorMatrixTree<T>& hmat);

typedef FactorMatrixTree<complex<double>> FMatrixTreecd;

typedef FactorMatrixTree<double> FMatrixTreed;
