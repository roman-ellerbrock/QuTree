//
// Created by Roman Ellerbrock on 3/13/20.
//

#ifndef CDVR_H
#define CDVR_H
#include "DVR/TDDVR.h"
#include "TreeOperators/Potential.h"
#include "DVR/DeltaVTree.h"
#include "DVR/cdvr_functions.h"
#include "TreeClasses/TensorTreeFunctions.h"
#include "TreeClasses/SymTensorTree.h"

class CDVR {
public:
	CDVR(const Tree& tree);
//	CDVR(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree, size_t part = 0);
	~CDVR() = default;

	void Update(const Wavefunction& Psi, const PotentialOperator& V, const Tree& tree,
		size_t part = 0, bool out = false, ostream& os = cout);

/*	void Update2(Wavefunction Psi, const PotentialOperator& V,
		const Tree& smalltree, size_t part = 0,
		bool out = false, ostream& os = std::cout);*/

//	Tensorcd Apply(Tensorcd Phi, const Matrixcd& invsq_rho, const Node& node) const;
	Tensorcd apply(Tensorcd Phi, const SpectralDecompositioncd& rho_x, const Node& node) const;
	Tensorcd applySym(Tensorcd Phi, const SpectralDecompositioncd& rho_x, const Node& node) const;

	void update(SymTensorTree& Psi, const PotentialOperator& V,
		const Tree& tree, size_t part = 0, bool out = false, ostream& os = std::cout);

	TDDVR tddvr_;

private:
	Tree ltree_;

	TensorTreecd Vnode_;
	MatrixTreed Vedge_;

	DeltaVTree deltaV_;

	TensorTreecd Cdown_;
	TensorTreecd Cup_;

	WorkMemorycd mem_;

};

#endif //CDVR_H
