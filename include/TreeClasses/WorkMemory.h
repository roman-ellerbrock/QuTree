//
// Created by Roman Ellerbrock on 6/16/21.
//

#ifndef WORKMEMORY_H
#define WORKMEMORY_H
#include "TensorTree.h"

template <typename T>
class WorkMemory {
public:
	WorkMemory() = default;
	~WorkMemory() = default;
	WorkMemory(const Tree& tree) {
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

	WorkMemory(const SparseTree& stree) {

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

	Tensor<T> work1_, work2_, work3_;
};

typedef WorkMemory<complex<double>> WorkMemorycd;
typedef WorkMemory<double> WorkMemoryd;

#endif //WORKMEMORY_H
