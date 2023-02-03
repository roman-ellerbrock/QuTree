//
// Created by Roman Ellerbrock on 1/5/23.
//

#ifndef BLOCKTREESHAPE_H
#define BLOCKTREESHAPE_H
#include "TreeClasses/Discrete/LabelTree.h"
#include "TreeClasses/Discrete/U1Symmetry.h"

using BlockTensorShape = vector<TensorShape>;
using BlockTreeShape = NodeAttribute<BlockTensorShape>;
using Partition = vector<size_t>;
using Partitions = vector<Partition>;

vector<const ConfigurationTensor<>*> mask(const vector<const BlockTensor*>& As, const Partition& partition) {
	vector<const ConfigurationTensor<>*> S;
	for (size_t l = 0; l < As.size(); ++l) {
		const BlockTensor& A = *As[l];
		size_t p = partition[l];
		/// if p doesn't exist, return an empty vector
		if (A.find(p) == A.end()) {
			S.clear();
			return S;
		}
 		S.push_back(&A.at(p));
	}
	return S;
}


class BlockTree {
public:
	BlockTree() = default;
	~BlockTree() = default;

	explicit BlockTree(const Tree& tree, const Range& range,
		const vector<Labels>& leaf_labels, size_t max_dim) {
		initialize(tree, range, leaf_labels, max_dim);
	}

	void initialize(const Tree& tree, const Range& range,
		const vector<Labels>& leaf_labels, size_t max_dim) {
		labels_.initialize(tree, leaf_labels, range);
		dim_.initialize(labels_, max_dim, tree);
	}

	BlockTensor bottomBlockTensor(const Labels& labels) {
		BlockTensor A;
		size_t i = 0;
		for (const Label& label: labels) {
			ConfigurationTensor<> c = {{i++}};
			A[label] = c;
		}
		return A;
	}

	ConfigurationTensor<> upperTensorUp(const Label& label, const Node& node) {
		ConfigurationTensor<> a;


		return a;
	}

	/// 1:
	/// 01, 10, | 00
	/// 00 | 01, 10
	/// 2:
	/// 11 | 00 (2+0)
	/// 00 | 11 (0+2)
	/// 01, 10 | 01, 10 (1+1)
	size_t nNeighbors(const Node& node, int hole) {
		size_t nNeigh = 0;
		for (size_t k = 0; k < node.nChildren(); ++k) {
			if (k == hole) { continue; }
			nNeigh++;
		}
		if (!node.isToplayer() && (node.parentIdx() != hole)) {
			nNeigh++;
		}
		return nNeigh;
	}

	BlockTensor randomBlockTensor(const Node& node, int hole,
		const Range& range, const Tree& tree) {
		const auto& neighbors = tree.neighbors(node, hole);
		auto numbers = labels_.up_.mask(neighbors);
		const Labels& target_labels = labels_.up_[node];
		const vector<const BlockTensor*> tensors = up_.mask(neighbors);
		BlockTensor A;
		for (const Label& label : target_labels) {
			/// decompose symmetry label into realizations
			const auto ps = partitions(numbers, label);
			/// for each partition
			for (const Partition p : ps) {
				/// get Configurationtensor for each Label
				auto blocks = mask(tensors, p);
//				cartesianProduct(A, tensors);
			}
		}
		/**
		 * for label in Labels:
		 *   ps = partitions(label, node, hole)
		 */
		return A;
	}

	void print(const Tree& tree) {
		labels_.print(tree);
	}

	LabelTree labels_;
	LabelDimensionTree dim_;
	NodeAttribute<BlockTensor> up_;
	NodeAttribute<BlockTensor> down_;
};

#endif //BLOCKTREESHAPE_H

/**
 * stuff todo:
 * - == for node
 * - rename SparseTree to TreeMask
 * - use SparseTree for dense as well
 * - every NodeAttribute has a SparseTree built in
 * 		- enables operator<<
 * - replace Matrix by Tensor
 * -
 *
*/

