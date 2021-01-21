//
// Created by Roman Ellerbrock on 4/3/20.
//

#ifndef MATRIXTENSORTREE_H
#define MATRIXTENSORTREE_H
#include "Core/Wavefunction.h"
#include "TreeClasses/MatrixTreeFunctions.h"

class MatrixTensorTree : public pair<Wavefunction, MatrixTreecd> {
public:
	MatrixTensorTree() = default;
	MatrixTensorTree(const Wavefunction& Psi, const Tree& tree, bool orthogonal);
	/**
	 * \brief All-normalized wavefunction representation with A\tilde's and B_inv's.
	 */

	void Initialize(const Wavefunction& Psi,
		const Tree& tree, bool orthogonal);

	const TensorTreecd& nodes()const { return first; }
	TensorTreecd& nodes() { return first; }

	const MatrixTreecd& edges()const { return second; }
	MatrixTreecd& edges() { return second; }

	const TensorTreecd& DensityWeighted() const { return nodes(); }

	void buildNodes(TensorTreecd Psi, const Tree& tree);

	void buildEdges(const Tree& tree);

	void buildFromWeighted(const Tree& tree);

	TensorTreecd TopDownNormalized(const Tree& tree) const;

	TensorTreecd BottomUpNormalized(const Tree& tree) const;

};

bool IsWorking_bottomup(const MatrixTensorTree& Psi, const Tree& tree, double eps);
bool IsWorking_topdown(const MatrixTensorTree& Psi, const Tree& tree, double eps);
bool IsWorking(const MatrixTensorTree& Psi, const Tree& tree, double eps = 1e-7);


#endif //MATRIXTENSORTREE_H
