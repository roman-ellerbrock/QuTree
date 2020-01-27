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
class HHoleMatrices: public SparseTreeStructuredObject<Matrixcd> {
public:
	HHoleMatrices(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrixcd>(active_, basis) {}

	HHoleMatrices(const MPO<T>& M, const TTBasis& basis)
		: SparseTreeStructuredObject<Matrixcd>(M, basis) {}

	~HHoleMatrices() = default;

	// Calculate Hole-Matrices
	void Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const HMatrices<T>& hmat, const TTBasis& basis);

	void Calculate(const TensorTree<T>& Psi,
		const HMatrices<T>& hmat, const TTBasis& basis) {
		Calculate(Psi, Psi, hmat, basis);
	}

	Tensorcd Apply(const Tensorcd& Phi, const Node& node) const;

	void print(TTBasis& basis, ostream& os = cout);
};


