//
// Created by Roman Ellerbrock on 2020-01-21.
//
#include "TensorTreeBasis.h"

TensorTreeBasis::TensorTreeBasis(const TensorTreeBasis& T)
	: tree(T.tree){
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

void TensorTreeBasis::Read(istream& file) {
	cout << "==================================================" << endl;
	cout << "=====          Basis Initialization           ====" << endl;
	cout << "==================================================" << endl << endl;

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

	cout << "==================================================" << endl;
	cout << "=====   - Basis initialized.                  ====" << endl;
	cout << "==================================================" << endl << endl;
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
	os << "List of Physical Coordinates:" << endl;
	for (size_t i = 0; i < this->nLeaves(); i++) {
		const Leaf& node = GetLeaf(i);
		node.info(os);
		os << endl;
	}
	os << endl;

	// ... and now for every logical node
	cout << "List of logical nodes:" << endl;
	cout << "nTotalNodes = " << nNodes() << endl;
	for (size_t i = 0; i < nNodes(); i++) {
		const Node& node = GetNode(i);
		node.info();
		node.TDim().print(cout);
		if (!node.IsToplayer()) {
			const Node& parent = node.Up();
			const TensorDim& parentdim = parent.TDim();
			parentdim.print(cout);
		}
		cout << endl;
	}
	cout << endl;
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
	cout << "Linearized " << counter + 1 << " nodes.\n";
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

