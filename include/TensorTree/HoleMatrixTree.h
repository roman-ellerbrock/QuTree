/**
 * \class HHoleMatrices
 *
 * \ingroup MCTDH-Matrices
 *
 * \brief This class represents the H-Hole-Matrices.
 *
 * The H-Hole-Matrices are objects that occur in the
 * equations of motion (EOM) of the MCTDH approach. A H-Hole-Matrix
 * is the mean-field operator representation of a MultiParticleOperator
 * for a given mctdhWavefunction.
 * */
#pragma once
#include "TreeStructuredObject.h"
#include "TensorTree.h"
#include "SingleParticleOperator.h"
#include "FactorMatrixTree.h"
#include "SparseTreeStructuredObject.h"

template<typename T>
class HoleMatrixTree: public SparseTreeStructuredObject<Matrix<T>> {
public:
	using SparseTreeStructuredObject<Matrix<T>>::Active;
	using SparseTreeStructuredObject<Matrix<T>>::operator[];
	using SparseTreeStructuredObject<Matrix<T>>::Initialize;
	using SparseTreeStructuredObject<Matrix<T>>::attributes;

	HoleMatrixTree(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(active_, basis) {}

	HoleMatrixTree(const MPO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(cast_to_vector_size_t(M.Modes()), basis) {
		Initialize(basis);
	}

	HoleMatrixTree(const TensorTree<T>& Psi, const FactorMatrixTree<T>& hmat,
		const MPO<T>& M, const TTBasis& basis);

	~HoleMatrixTree() = default;

	void Initialize(const TTBasis& basis) override;

	// Calculate Hole-Matrices
	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const FactorMatrixTree<T>& hmat, const TTBasis& basis);

	void Calculate(const TensorTree<T>& Psi,
		const FactorMatrixTree<T>& hmat, const TTBasis& basis) {
		Calculate(Psi, Psi, hmat, basis);
	}

	Tensor<T> Apply(const Tensor<T>& Phi, const Node& node) const;

	/// I/O
	void print(TTBasis& basis, ostream& os = cout);
	void Write(ostream& os) const;
	void Write(const string& filename) const;
	void Read(istream& is);
	void Read(const string& filename);
};

template<typename T>
ostream& operator>>(ostream& os, const HoleMatrixTree<T>& hmat);

template<typename T>
istream& operator<<(istream& is, HoleMatrixTree<T>& hmat);

typedef HoleMatrixTree<complex<double>> HHoleMatricescd;

typedef HoleMatrixTree<double> HHoleMatricesd;



