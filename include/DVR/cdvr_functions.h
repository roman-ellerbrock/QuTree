//
// Created by Roman Ellerbrock on 3/18/20.
//

#ifndef CDVR_FUNCTIONS_H
#define CDVR_FUNCTIONS_H
#include "TreeClasses/MatrixTree.h"
#include "DVR/DeltaVTree.h"
#include "DVR/MatrixTensorTree.h"
#include "DVR/TreeGrids.h"

namespace cdvr_functions {

	void fillXNode(Vectord& X, vector<size_t> idx, const TreeGrids& grids, const TreeGrids& holegrids,
		const Node& node);

	void fillXEdge(Vectord& X, vector<size_t> idx, const TreeGrids& grids, const TreeGrids& holegrids,
		const Node& node);

	void calculateDeltaVs(DeltaVTree& deltaVs,
		const TensorTreecd& Cup, const TensorTreecd& Cdown,
		const TensorTreecd& Vnodes, const MatrixTreed& Vedges,
		const Tree& tree);

//	TensorTreecd Apply(const Wavefunction& Psio, const TensorTreecd& Cdown,
//		const TensorTreecd& V, const DeltaVTree& DeltaVs, const Tree& tree);

	void apply(Tensorcd& VPhi, const Tensorcd& Xi, const Tensorcd& V, const TensorTreecd& Cdown,
		const DeltaVTree& deltaVs, const Node& node, const WorkMemorycd& mem);
}

#endif //CDVR_FUNCTIONS_H
