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
#include "HMatrices.h"
#include "SparseTreeStructuredObject.h"

template<typename T>
class HHoleMatrices: public SparseTreeStructuredObject<Matrix<T>> {
public:
	using SparseTreeStructuredObject<Matrix<T>>::Active;
	using SparseTreeStructuredObject<Matrix<T>>::operator[];
	using SparseTreeStructuredObject<Matrix<T>>::Initialize;
	using SparseTreeStructuredObject<Matrix<T>>::attributes;

	HHoleMatrices(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(active_, basis) {}

	HHoleMatrices(const MPO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrix<T>>(cast_to_vector_size_t(M.Modes()), basis) {
		Initialize(basis);
	}

	HHoleMatrices(const TensorTree<T>& Psi, const HMatrices<T>& hmat,
		const MPO<T>& M, const TTBasis& basis);

	~HHoleMatrices() = default;

	void Initialize(const TTBasis& basis) override;

	// Calculate Hole-Matrices
	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const HMatrices<T>& hmat, const TTBasis& basis);

	void Calculate(const TensorTree<T>& Psi,
		const HMatrices<T>& hmat, const TTBasis& basis) {
		Calculate(Psi, Psi, hmat, basis);
	}

	Tensor<T> Apply(const Tensor<T>& Phi, const Node& node) const;

	void print(TTBasis& basis, ostream& os = cout);
};

typedef HHoleMatrices<complex<double>> HHoleMatricescd;
typedef HHoleMatrices<double> HHoleMatricesd;



