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
#include "Matrix.h"
#include "MultiParticleOperator.h"
#include "TensorTreeBasis.h"
#include "TensorTree.h"
#include "SparseTreeStructuredObject.h"
#include "FactorMatrix.h"

vector<size_t> cast_to_vector_size_t(const vector<int>& a);

template<typename T>
class HMatrices: public SparseTreeStructuredObject<FactorMatrix<T>> {
public:
	using SparseTreeStructuredObject<FactorMatrix<T>>::Active;
	using SparseTreeStructuredObject<FactorMatrix<T>>::operator[];
	using SparseTreeStructuredObject<FactorMatrix<T>>::Initialize;
	using SparseTreeStructuredObject<FactorMatrix<T>>::attributes;

	HMatrices(const MPO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<FactorMatrix<T>>(cast_to_vector_size_t(M.Modes()), basis) {
		Initialize(basis);
	}

	HMatrices(const TensorTree<T>& Psi, const MPO<T>& M, const TTBasis& basis)
		: HMatrices(M, basis) {
		Calculate(Psi, M, basis);
	}

	HMatrices(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<FactorMatrix<T>>(active_, basis) {
		Initialize(basis);
	}

	HMatrices(const MPO<T>& M, const TTBasis& basis, const string& filename)
		: HMatrices(M, basis) {
		Read(filename);
	}

	~HMatrices() = default;

	void Initialize(const TTBasis& basis) override;

	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const MPO<T>& M, const TTBasis& basis);

	void Calculate(const TensorTree<T>& Psi,
		const MPO<T>& M, const TTBasis& basis) {
		Calculate(Psi, Psi, M, basis);
	}

	void CalculateLayer(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const MPO<T>& M, const Node& node);

	Tensor<T> Apply(const Tensor<T>& Phi, const MPO<T>& M, const Node& node) const;

	Tensor<T> ApplyUpper(Tensor<T> Phi, const Node& node) const;

	Tensor<T> ApplyHole(Tensor<T> Phi, const Node& hole_node) const;

	/// I/O
	void print(const TTBasis& basis, ostream& os = cout) const;
	void Write(ostream& os) const;
	void Write(const string& filename) const;
	void Read(istream& is);
	void Read(const string& filename);

protected:
	FactorMatrix<T> CalculateUpper(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const Node& node);

	FactorMatrix<T> CalculateBottom(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const MPO<T>& M, const Node& node, const Leaf& phys);
};

template <typename T>
ostream& operator>>(ostream& os, const HMatrices<T>& hmat);

template <typename T>
istream& operator<<(istream& is, HMatrices<T>& hmat);

typedef HMatrices<complex<double>> HMatricescd;

typedef HMatrices<double> HMatricesd;
