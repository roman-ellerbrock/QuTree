//
// Created by Roman Ellerbrock on 2/2/20.
//

#include "TreeShape/TreeFactory.h"

namespace TreeFactory {

	Node InitTrain(const Node& bottom, size_t dimNodes) {
		Node train;
		vector<size_t> dims;
		train.push_back(bottom);
		dims.push_back(dimNodes);
		train.push_back(bottom);
		dims.push_back(dimNodes);
		dims.push_back(dimNodes);
		TensorShape tdim(dims);
		train.shape() = tdim;
		return train;
	}

	Node TrainLayer(const Node& old_train, const Node& bottom, size_t dimNodes) {
		Node train;
		train.push_back(bottom);
		train.push_back(old_train);
		TensorShape tdim({dimNodes, dimNodes, dimNodes});
		train.shape() = tdim;
		return train;
	}

	Tree UnbalancedTree(size_t nLeaves, size_t dimLeaves, size_t dimNodes, size_t leafType) {
		size_t mode = 0;
		size_t leafSubtype = 0;
		PhysPar par;
		Leaf leaf(dimLeaves, mode, leafType, leafSubtype, par);
		Node bottom(leaf, dimNodes);
		Node train = InitTrain(bottom, dimNodes);

		for (size_t k = 0; k < nLeaves - 2; ++k) {
			train = TrainLayer(train, bottom, dimNodes);
		}

		train.setParent(nullptr);
		auto& tdim_ = train.shape();
		tdim_.setDimension(1, tdim_.lastIdx());
		Tree tree;
		tree.SetRoot(train);
		tree.ResetLeafModes();
		return tree;
	}

	vector<Node> Partition(const vector<Node>& nodes,
		size_t n_partition, size_t dim_node) {
		/// This is a helper function to create close-to balanced trees.
		/// It adds a layer to a vector of nodes
		vector<Node> groups;
		size_t n_loop = nodes.size() / n_partition;
		for (size_t k = 0; k < n_loop; ++k) {
			Node p;
			vector<size_t> dims;
			for (size_t l = 0; l < n_partition; ++l) {
				p.push_back(nodes[k * n_partition + l]);
				dims.push_back(dim_node);
			}
			dims.push_back(dim_node);
			TensorShape tensordim(dims);
			p.shape() = tensordim;
			groups.emplace_back(p);
		}
		size_t n_rest = nodes.size() % n_partition;
		for (size_t r = 0; r < n_rest; ++r) {
			groups.push_back(nodes[n_loop * n_partition + r]);
		}
		return groups;
	}

	Tree BalancedTree(size_t num_leaves,
		size_t dim_leaves, size_t dim_nodes) {
		/// Create close-to-balanced Tree
		size_t leaf_type = 6;
		size_t mode = 0;
		size_t leaf_subtype = 0;
		PhysPar par;
		Leaf leaf(dim_leaves, mode, leaf_type, leaf_subtype, par);

//		size_t dim_now = min(dim_nodes, dim_leaves);
		size_t dim_now = dim_nodes;
		Node bottom(leaf, dim_now);
		vector<Node> nodes;
		for (size_t k = 0; k < num_leaves; ++k) {
			nodes.push_back(bottom);
		}
		size_t count = 0;
		while (nodes.size() > 1) {
			nodes = Partition(nodes, 2, dim_nodes);
			count++;
			if (count > 100) {
				cerr << "Error while partitioning TensorTreeBasis in constructor.\n";
				exit(1);
			}
		}
		auto& tdim = nodes.front().shape();
		tdim.setDimension(1, tdim.lastIdx());
		Tree tree;
		tree.SetRoot(nodes.front());
		tree.ResetLeafModes();
		return tree;
	}
}


