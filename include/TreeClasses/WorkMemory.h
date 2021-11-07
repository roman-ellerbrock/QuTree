//
// Created by Roman Ellerbrock on 6/16/21.
//

#ifndef WORKMEMORY_H
#define WORKMEMORY_H
#include "TensorTree.h"

template<typename T>
class WorkMemory {
public:
	WorkMemory() = default;
	~WorkMemory() = default;

	explicit WorkMemory(const Tree& tree) {
		TensorShape shape({1});
		for (const Node& node : tree) {
			if (node.shape().totalDimension() > shape.totalDimension()) {
				shape = node.shape();
			}
		}
		work1_ = Tensor<T>(shape);
		work2_ = Tensor<T>(shape);
		work3_ = Tensor<T>(shape);
	}

	explicit WorkMemory(const SparseTree& stree) {

		int sub_topnode = stree.size() - 1;
		/// Find largest TensorShape in tree and allocate work memory
		TensorShape shape({1});
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = stree.node(n);
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				if (parent.shape().totalDimension() > shape.totalDimension()) {
					shape = parent.shape();
				}
			}
		}
		work1_ = Tensor<T>(shape);
		work2_ = Tensor<T>(shape);
		work3_ = Tensor<T>(shape);
	}

	explicit WorkMemory(const Node& node) {
		TensorShape shape = node.shape();
		if (!node.isToplayer()) {
			const Node& parent = node.parent();
			if (parent.shape().totalDimension() > shape.totalDimension()) {
				shape = parent.shape();
			}
		}
		work1_ = Tensor<T>(shape);
		work2_ = Tensor<T>(shape);
		work3_ = Tensor<T>(shape);
	}

	mutable Tensor<T> work1_; /// state of workmem is irrelevant,
	mutable Tensor<T> work2_; /// marked mutable to allow const even if work mem is member of classes
	mutable Tensor<T> work3_;
};

typedef WorkMemory<complex<double>> WorkMemorycd;

typedef WorkMemory<double> WorkMemoryd;

#endif //WORKMEMORY_H
