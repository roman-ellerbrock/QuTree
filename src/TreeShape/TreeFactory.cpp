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

	Tree unbalancedTree(size_t nLeaves, size_t dimLeaves, size_t dimNodes, size_t leafType) {
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
		tree.setRoot(train);
		tree.resetLeafModes();
		return tree;
	}

	vector<Node> partition(const vector<Node>& nodes,
		size_t n_partition, size_t dim_node) {
		/// This is a helper function to create close-to balanced trees.
		/// It adds a layer to a vector of nodes
		vector<Node> groups;
		size_t n_loop = nodes.size() / n_partition;
		for (size_t k = 0; k < n_loop; ++k) {
			Node p;
			vector<size_t> dims;
			size_t product = 1;
			for (size_t l = 0; l < n_partition; ++l) {
				const Node& child = nodes[k * n_partition + l];
				p.push_back(child);
				size_t dim_now = child.shape().lastDimension();
				product *= dim_now;
				dims.push_back(dim_now);
			}
			size_t dim_parent = min(product, dim_node);
			dims.push_back(dim_parent);
			TensorShape tensordim(dims);
			p.shape() = tensordim;
			groups.emplace_back(p);
		}

		size_t n_rest = nodes.size() % n_partition;
		if ((n_loop == 0) && (n_rest > 0)) {
			/// if nothing would be done, wrap it up
			Node p;
			vector<size_t> dims;
			size_t product = 1;
			for (size_t l = 0; l < n_rest; ++l) {
				const Node& child = nodes[l];
				p.push_back(child);
				size_t dim_now = child.shape().lastDimension();
				product *= dim_now;
				dims.push_back(dim_now);
			}
			size_t dim_parent = min(product, dim_node);
			dims.push_back(dim_parent);
			TensorShape tensordim(dims);
			p.shape() = tensordim;
			groups.emplace_back(p);
		} else {
			/// Just append and deal with dangling nodes later
			for (size_t r = 0; r < n_rest; ++r) {
				groups.push_back(nodes[n_loop * n_partition + r]);
			}
		}

		return groups;
	}

	vector<Node> bottomlayerNodes(size_t num_leaves, size_t dim_leaves, size_t dim_nodes) {
		/// Hardcoded leaf-parameters
		size_t leaf_type = 6;
		size_t mode = 0;
		size_t leaf_subtype = 0;
		PhysPar par;
		Leaf leaf(dim_leaves, mode, leaf_type, leaf_subtype, par);

		/// Creat bottomlayer nodes manually
		size_t dim_now = min(dim_nodes, dim_leaves);
		Node bottom(leaf, dim_now);
		vector<Node> nodes;
		for (size_t k = 0; k < num_leaves; ++k) {
			nodes.push_back(bottom);
		}
		return nodes;
	}

	Tree expandNodes(Tree tree) {
		for (auto it = tree.begin(); it != tree.end(); ++it) {
			Node& node = *it;
			const TensorShape& shape = node.shape();
			if ((!node.isBottomlayer()) && (!node.isToplayer())) {
				if (shape.lastDimension() == shape.lastBefore()) {
					/// expand
					Node& parent = node.parent();
					parent.expandChild(node.childIdx());
					it--;
					tree.update();
				}
			}
		}
		return tree;
	}

	Tree balancedTree(size_t num_leaves,
		size_t dim_leaves, size_t dim_nodes) {
		/**
		 * \brief This functions creates a close-to-balanced Tree
		 * \@param num_leaves number of leaves in the tree
		 * \@param dim_leaves dimension of basis at leaves (primitive or physical basis)
		 * \@param dim_nodes dimnesion at higher nodes, "virtual bond dimension" or "number of SPFs"
		 * \@return the generated (close-to-)balanced tree
		 */

		/// Cover leaves in bottomlayer nodes
		auto nodes = bottomlayerNodes(num_leaves, dim_leaves, dim_nodes);

		/// Add layer after layer until only one node is left
		size_t count = 0;
		while (nodes.size() > 1) {
			nodes = partition(nodes, 2, dim_nodes);
			count++;
			if (count > 100) {
				cerr << "Error while partitioning TensorTreeBasis in constructor.\n";
				exit(1);
			}
		}
		/// Set number of wavefunctions to 1 and asign root-node
		auto& tdim = nodes.front().shape();
		tdim.setDimension(1, tdim.lastIdx());
		Tree tree;
		tree.setRoot(nodes.front());
		tree.resetLeafModes();

		/// Expand nodes that perform no contraction
		tree = expandNodes(tree);

		return tree;
	}

	Tree operatorTree(const Tree& tree) {
		Tree otree(tree);
		for (Node& node : otree) {
			if (node.isBottomlayer()) {
				TensorShape& shape = node.shape();
				size_t dim = shape.lastBefore();
				size_t ldim = shape.lastDimension();
				TensorShape newshape({dim * dim, ldim * ldim});
				node.shape() = newshape;

				Node& parent = node.parent();
				TensorShape& pshape = parent.shape();
				pshape.setDimension(ldim * ldim, node.childIdx());
			}
		}
		return otree;
	}
}


