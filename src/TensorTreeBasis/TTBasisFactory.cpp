//
// Created by Roman Ellerbrock on 2/2/20.
//

#include "TreeHandling/TTBasisFactory.h"

namespace TTBasisFactory {

	Node InitTrain(const Node& bottom, size_t dimNodes) {
		Node train;
		vector<size_t> dims;
		train.push_back(bottom);
		dims.push_back(dimNodes);
		train.push_back(bottom);
		dims.push_back(dimNodes);
		TensorDim tdim(dims, dimNodes);
		train.TDim() = tdim;
		return train;
	}

	Node TrainLayer(const Node& old_train, const Node& bottom, size_t dimNodes) {
		Node train;
		train.push_back(bottom);
		train.push_back(old_train);
		TensorDim tdim({dimNodes, dimNodes}, dimNodes);
		train.TDim() = tdim;
		return train;
	}

	Tree TensorTrain(size_t nLeaves, size_t dimLeaves, size_t dimNodes, size_t leafType) {
		size_t mode = 0;
		size_t leafSubtype = 0;
		PhysPar par;
		Leaf leaf(dimLeaves, mode, leafType, leafSubtype, par);
		Node bottom(leaf, dimNodes);
		Node train = InitTrain(bottom, dimNodes);

		for (size_t k = 0; k < nLeaves - 2; ++k) {
			train = TrainLayer(train, bottom, dimNodes);
		}

		train.SetUp(nullptr);
		auto& tdim_ = train.TDim();
		tdim_.SetNumTensor(1);
		Tree basis;
		basis.SetRoot(train);
		basis.Update();
		ResetLeafModes(basis);
//		tree = move(train);
//		tree.SetUp(nullptr);
//		tree.UpdatePosition(NodePosition());
//		Update();
//		ResetLeafModes(*this);
//		basis.ResetLeafModes(basis);
		return basis;
	}
}


