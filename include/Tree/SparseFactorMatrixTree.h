#pragma once
#include "Core/Matrix.h"
#include "TensorTree.h"
#include "SparseTreeStructuredObject.h"
#include "MultiLeafOperator.h"
#include "Core/FactorMatrix.h"


vector<size_t> cast_to_vector_size_t(const vector<int>& a);

template<typename T>
class SparseFactorMatrixTree: public SparseTreeStructuredObject<FactorMatrix<T>>
/**
 * \class FactorMatrixTree
 *
 * \ingroup Tree
 *
 * \brief This class represents the FactorMatrices for Trees.
 *
 * The FactorMatrices are the result of Tensor-products of tensor trees.
 * In a physical context, the hole-matrices are representation of
 * mean-field operators when working with tensor tree wavefunctions.
 * */
	{
public:
	using SparseTreeStructuredObject<FactorMatrix<T>>::Active;
	using SparseTreeStructuredObject<FactorMatrix<T>>::operator[];
	using SparseTreeStructuredObject<FactorMatrix<T>>::Initialize;
	using SparseTreeStructuredObject<FactorMatrix<T>>::attributes_;

	/// Create FactorMatrices for relevant nodes when representing an operator
	SparseFactorMatrixTree(const MLO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<FactorMatrix<T>>(cast_to_vector_size_t(M.Modes()), basis) {
		Initialize(basis);
	}

	/// Create and calculate FactorMatrixTree for an operator
	SparseFactorMatrixTree(const TensorTree<T>& Psi, const MLO<T>& M, const TTBasis& basis)
		: SparseFactorMatrixTree(M, basis) {
		Calculate(Psi, M, basis);
	}

	/// Create FactorMatrices for externally marked nodes
	SparseFactorMatrixTree(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<FactorMatrix<T>>(active_, basis) {
		Initialize(basis);
	}

	/// Read FactorMatrixTree from file
	SparseFactorMatrixTree(const MLO<T>& M, const TTBasis& basis, const string& filename)
		: SparseFactorMatrixTree(M, basis) {
		Read(filename);
	}

	/// Default destructor
	~SparseFactorMatrixTree() = default;

	/// Fill FactorMatrices at internally marked nodes. TreeMarker must be initialized
	void Initialize(const TTBasis& basis) override;

	/// Calculate (cross-)TreeMatrix-representation of an operator
	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const MLO<T>& M, const TTBasis& basis);

	/// Calculate TreeMatrix-representation of an operator
	void Calculate(const TensorTree<T>& Psi,
		const MLO<T>& M, const TTBasis& basis) {
		Calculate(Psi, Psi, M, basis);
	}

	/// Calculate matrices locally at a node
	void CalculateLayer(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const MLO<T>& M, const Node& node);

	/// Apply factor matrices locally
	Tensor<T> Apply(const Tensor<T>& Phi, const MLO<T>& M, const Node& node) const;

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
		const MLO<T>& M, const Node& node, const Leaf& phys);
};

template<typename T>
ostream& operator>>(ostream& os, const SparseFactorMatrixTree<T>& hmat);

template<typename T>
istream& operator<<(istream& is, SparseFactorMatrixTree<T>& hmat);

typedef SparseFactorMatrixTree<complex<double>> SparseFactorMatrixTreecd;

typedef SparseFactorMatrixTree<double> SparseFactorMatrixTreed;
