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
#include "SingleParticleOperator.h"
#include "MultiParticleOperator.h"
#include "TensorTreeBasis.h"
#include "TensorTree.h"
#include "SingleParticleOperator.h"
#include "SparseTreeStructuredObject.h"

template<typename T>
class HMatrices: public SparseTreeStructuredObject<SPO<T>> {
public:
	using SparseTreeStructuredObject<SPO<T>>::Active;
	using SparseTreeStructuredObject<SPO<T>>::operator[];
/*	HMatrices(const MPO<T>& M,
		const TTBasis& basis)
		: SparseTreeStructuredObject<SPO<T>>(M, basis) {}
*/
	HMatrices(shared_ptr<TreeMarker>& active_, const TTBasis& basis)
		: SparseTreeStructuredObject<SPO<T>>(active_, basis) {}

	~HMatrices() = default;

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

	void print(const TTBasis& basis, ostream& os = cout) const;

protected:
	SPO<T> CalculateUpper(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const Node& node);

	SPO<T> CalculateBottom(const Tensor<T>& Bra, const Tensor<T>& Ket,
		const MPO<T>& M, const Node& node, const Leaf& phys);
};


