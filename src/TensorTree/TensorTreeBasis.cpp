//
// Created by Roman Ellerbrock on 2020-01-21.
//
#include "TensorTreeBasis.h"

TensorTreeBasis::TensorTreeBasis(const TensorTreeBasis& T)
	: tree(T.tree) {
	Update();
}

TensorTreeBasis::TensorTreeBasis(TensorTreeBasis&& T) noexcept {
	tree = move(T.tree);
	linearizedNodes_ = move(T.linearizedNodes_);
	linearizedLeaves_ = move(T.linearizedLeaves_);
}

TensorTreeBasis& TensorTreeBasis::operator=(const TensorTreeBasis& T) {
	*this = TTBasis(T);
	return *this;
}

TensorTreeBasis& TensorTreeBasis::operator=(TensorTreeBasis&& T) noexcept {
	std::swap(tree, T.tree);
	std::swap(linearizedNodes_, T.linearizedNodes_);
	std::swap(linearizedLeaves_, T.linearizedLeaves_);
	return *this;
}

TensorTreeBasis::TensorTreeBasis(const string& filename) {
	Read(filename);
}

TensorTreeBasis::TensorTreeBasis(istream& is) {
	Read(is);
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
		TensorDim tensordim(dims, dim_node);
		p.TDim() = tensordim;
		groups.emplace_back(p);
	}
	size_t n_rest = nodes.size() % n_partition;
//	auto& last = groups.back();
	for (size_t r = 0; r < n_rest; ++r) {
//		last.push_back(nodes[n_loop * n_partition + r]);
		groups.push_back(nodes[n_loop * n_partition + r]);
	}
	return groups;
}

void ResetLeafModes(TensorTreeBasis& basis) {
	size_t n_modes = basis.nLeaves();
	assert(n_modes > 0);
	size_t mode = n_modes - 1;
	for (Node& node : basis) {
		if (node.IsBottomlayer()) {
			Leaf& leaf = node.PhysCoord();
			leaf.Mode() = mode--;
		}
	}
}

TensorTreeBasis::TensorTreeBasis(size_t order,
	/// Create close-to-balanced Tree
	size_t dim_leaves, size_t dim_nodes) {
	size_t leaf_type = 6;
	size_t mode = 0;
	size_t leaf_subtype = 0;
	PhysPar par;
	Leaf leaf(dim_leaves, mode, leaf_type, leaf_subtype, par);

	Node bottom(leaf, dim_nodes);
	vector<Node> nodes;
	for (size_t k = 0; k < order; ++k) {
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
	tree = move(nodes.front());
	tree.SetUp(nullptr);
	auto& tdim = tree.TDim();
	tdim.setntensor(1);
	tree.UpdatePosition(NodePosition());
	Update();
	ResetLeafModes(*this);
}

void TensorTreeBasis::Read(istream& file) {
	// feed linearizedLeaves_ and logical block with references
	tree.Initialize(file, nullptr, NodePosition());
	Update();

	// Add new PhysPar for every physical coordinate
	for (int i = 0; i < linearizedLeaves_.size(); i++) {
		// Set parameters
		PhysPar par(file);
		linearizedLeaves_[i].SetPar(par);

		// Initialize primitive grid (HO, FFT, Legendre, ...)
		PrimitiveBasis& primitivebasis = linearizedLeaves_[i].PrimitiveGrid();
		primitivebasis.Initialize(par.Omega(), par.R0(), par.WFR0(), par.WFOmega());
	}
}

void TensorTreeBasis::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

void TensorTreeBasis::Write(ostream& os) const {
	tree.Write(os);
	linearizedLeaves_.Write(os);
}

void TensorTreeBasis::info(ostream& os) const {
	os << "List of Leaves:" << endl;
	for (size_t i = 0; i < this->nLeaves(); i++) {
		const Leaf& node = GetLeaf(i);
		node.info(os);
		os << endl;
	}
	os << endl;

	// ... and now for every logical node
	os << "List of upper nodes:" << endl;
	for (int i = nNodes() - 1; i >= 0; i--){
		const Node& node = GetNode(i);
		node.info();
		node.TDim().print(os);
		os << endl;
	}
	os << "Number of Nodes = " << nNodes() << endl;
}

Leaf& TensorTreeBasis::GetLeaf(size_t i) {
	return linearizedLeaves_[i];
}

const Leaf& TensorTreeBasis::GetLeaf(size_t i) const {
	return linearizedLeaves_[i];
}

Node& TensorTreeBasis::GetNode(size_t i) {
	return linearizedNodes_[i];
}

const Node& TensorTreeBasis::GetNode(size_t i) const {
	return linearizedNodes_[i];
}

void TensorTreeBasis::ExpandNode(Node& node) {
	assert(!node.IsToplayer());
	assert(!node.IsBottomlayer());

	Node& parent = node.Up();
	size_t childIdx = node.ChildIdx();
	parent.ExpandChild(childIdx);
	LinearizeNodes();
}

void TensorTreeBasis::Update() {
	// Tree is assumed to be updated, but the rest not:
	// Update everything
	tree.Update(NodePosition());
	LinearizeLeaves();
	LinearizeNodes();
}

void TensorTreeBasis::ReplaceNode(Node& old_node, Node& new_node) {
	// The old node must not be the toplayer node, otherwise change the
	// whole tree
	assert(!old_node.IsToplayer());

	// Replace the node
	Node& parent = old_node.Up();
	parent.Replace(new_node, old_node.ChildIdx());

	Node& topnode = TopNode();
	topnode.UpdatePosition(NodePosition());
	topnode.Updatennodes();
	LinearizeLeaves();
	LinearizeNodes();
}

void TensorTreeBasis::LinearizeNodes() {
	// block has to be cleared, because logical block
	// must be resistant to re-feed (important for e.g. expand node)
	// This routine adds every mctdh node to the logical block
	linearizedNodes_.clear();
	int counter = 0;
	for (int i = 0; i < nTotalNodes(); i++) {
		AbstractNode& abstract_node = nextNode();
		if (abstract_node.NodeType() == 1) {
			auto& node = (Node&) (abstract_node);
			node.SetAddress(counter);
			counter++;
			linearizedNodes_.push_back(node);
		}
	}
}

void TensorTreeBasis::LinearizeLeaves() {
	// This routine attends physical coordinates to the linearizedLeaves_ block
	linearizedLeaves_.clear();
	linearizedLeaves_.resizeaddress(nLeaves());

	for (int i = 0; i < nTotalNodes(); i++) {
		AbstractNode& abstract_node = nextNode();
		// If this node is a physical mode push it back
		if (abstract_node.NodeType() == 0) {
			auto& leaf = (Leaf&) (abstract_node);
			linearizedLeaves_.push_back(leaf);
			// @TODO: Check leaf-index mapping
//			reference_wrapper<Leaf> newphysmode(physnode);
//			linearizedLeaves_(physnode.Mode()) = newphysmode;
		}
	}
}

