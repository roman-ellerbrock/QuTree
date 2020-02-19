//
// Created by Roman Ellerbrock on 2020-01-21.
//
#include "TreeShape/Tree.h"

Tree::Tree(const Tree& T)
	: root_(T.root_) {
	Update();
}

Tree::Tree(Tree&& T) noexcept {
	root_ = move(T.root_);
	Update();
}

Tree& Tree::operator=(const Tree& T) {
	*this = Tree(T);
	return *this;
}

Tree& Tree::operator=(Tree&& T) noexcept {
	root_ = move(T.root_);
	Update();
	return *this;
}

Tree::Tree(const string& filename) {
	Read(filename);
}

Tree::Tree(istream& is) {
	Read(is);
}

void Tree::ResetLeafModes() {
	size_t n_modes = this->nLeaves();
	assert(n_modes > 0);
	int mode = n_modes - 1;
	for (Node& node : *this) {
		if (node.isBottomlayer()) {
			Leaf& leaf = node.getLeaf();
			leaf.Mode() = mode--;
		}
	}
	this->Update();
}

void Tree::ReindexLeafModes(map<size_t, size_t> Map) {
	for (Node& node : *this) {
		if (node.isBottomlayer()) {
			Leaf& leaf = node.getLeaf();
			leaf.Mode() = Map[leaf.Mode()];
		}
	}
}

void Tree::ExpandNode(Node& node) {
	assert(!node.isToplayer());
	assert(!node.isBottomlayer());

	Node& parent = node.parent();
	size_t childIdx = node.childIdx();
	parent.expandChild(childIdx);
	LinearizeNodes();
}

void Tree::Update() {
	// Tree is assumed to be updated, but the rest not:
	// Update everything
	root_.update(NodePosition());
	LinearizeLeaves();
	LinearizeNodes();
}

void Tree::ReplaceNode(Node& old_node, Node& new_node) {
	// The old node must not be the toplayer node, otherwise change the
	// whole tree
	assert(!old_node.isToplayer());

	// Replace the node
	Node& parent = old_node.parent();
	parent.Replace(new_node, old_node.childIdx());

	Node& topnode = TopNode();
	topnode.UpdatePosition(NodePosition());
	topnode.Updatennodes();
	LinearizeLeaves();
	LinearizeNodes();
}

void Tree::LinearizeNodes() {
	// block has to be cleared, because logical block
	// must be resistant to re-feed (important for e.g. expand node)
	// This routine adds every mctdh node to the logical block
	linearizedNodes_.clear();
	int counter = 0;
	for (int i = 0; i < nTotalNodes(); i++) {
		AbstractNode& abstract_node = nextNode();
		if (abstract_node.type() == 1) {
			auto& node = (Node&) (abstract_node);
			node.SetAddress(counter);
			counter++;
			linearizedNodes_.push_back(node);
		}
	}
}

void Tree::LinearizeLeaves() {
	// This routine attends physical coordinates to the linearizedLeaves_ block
	linearizedLeaves_.clear();
	linearizedLeaves_.resizeaddress(nLeaves());

	for (int i = 0; i < nTotalNodes(); i++) {
		AbstractNode& abstract_node = nextNode();
		// If this node is a physical mode push it back
		if (abstract_node.type() == 0) {
			auto& leaf = (Leaf&) (abstract_node);
			linearizedLeaves_.push_back(leaf);
			// @TODO: Check leaf-index mapping
//			reference_wrapper<Leaf> newphysmode(physnode);
//			linearizedLeaves_(physnode.Mode()) = newphysmode;
		}
	}
}

/// I/O

void Tree::Read(istream& file) {
	// feed linearizedLeaves_ and logical block with references
	root_.Initialize(file, nullptr, NodePosition());
	Update();

	// Add new PhysPar for every physical coordinate
	for (int i = 0; i < linearizedLeaves_.size(); i++) {
		// Set parameters
		PhysPar par(file);
		linearizedLeaves_[i].SetPar(par);

		// Initialize primitive grid (HO, FFT, Legendre, ...)
		LeafInterface& primitivebasis = linearizedLeaves_[i].PrimitiveGrid();
		primitivebasis.Initialize(par.Omega(), par.R0(), par.WFR0(), par.WFOmega());
	}
}

void Tree::Read(const string& filename) {
	ifstream is(filename);
	Read(is);
}

void Tree::Write(ostream& os) const {
	root_.Write(os);
}

ostream& operator<<(ostream& os, Tree& basis) {
	basis.Write(os);
	return os;
}

istream& operator<<(istream& is, Tree& basis) {
	basis.Read(is);
	return is;
}

void Tree::info(ostream& os) const {
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
		node.shape().print(os);
		os << endl;
	}
	os << "Number of Nodes = " << nNodes() << endl;
}

bool Tree::IsWorking() {

	bool works = true;
	int counter = 0;
	for (int i = 0; i < nTotalNodes(); i++) {
		const AbstractNode& abstract_node = nextNode();
		if (abstract_node.type() == 1) {
			auto& node = (Node&) (abstract_node);
			// 1.) Check global address
			// Do not break here to leave nodes in a valid state
			if (counter != node.Address()) { works = false; }
			counter++;

			if (!node.isBottomlayer()) {
				for (size_t k = 0; k < node.nChildren(); ++k) {
					const Node& child = node.child(k);
					if (&node != &child.parent()) {
						cerr << "Connectivity between child and parent is broken." << endl;
						return false;
					}
				}
			}
		}
	}
	if (!works) {
		cerr << "Corrupted address in tree." << endl;
		return false;
	}

	/// Check linearized Nodes
	if (nNodes() != linearizedNodes_.size()) {
		cerr << "linearizedNodes_ size does not match tree-size" << endl;
		return false;
	}
	counter = 0;
	for (int i = 0; i < nTotalNodes(); i++) {
		const AbstractNode& abstract_node = nextNode();
		if (abstract_node.type() == 1) {
			auto& node = (Node&) (abstract_node);
			// Do not break here to leave nodes in a valid state
			if (&node != &linearizedNodes_[counter].get()) { works = false; }
			counter++;
		}
	}
	if (!works) {
		cerr << "Corrupted linearizedNodes_. Missing Update()?" << endl;
		return false;
	}
	if (!linearizedNodes_.back().get().isToplayer()) {
		cerr << "Last node does not fulfill top-criterium." << endl;
		return false;
	}

	for (int i = 0; i < nTotalNodes(); i++) {
		AbstractNode& abstract_node = nextNode();
		// If this node is a physical mode push it back
		if (abstract_node.type() == 0) {
			auto& leaf = (Leaf&) (abstract_node);
			if (&leaf != &linearizedLeaves_[leaf.Mode()]) { works = false; }
		}
	}
	if (!works) {
		cerr << "Corrupted linearizedLeaves." << endl;
		return false;
	}

	return true;
}

Leaf& Tree::GetLeaf(size_t i) {
	return linearizedLeaves_[i];
}

const Leaf& Tree::GetLeaf(size_t i) const {
	return linearizedLeaves_[i];
}

Node& Tree::GetNode(size_t i) {
	return linearizedNodes_[i];
}

const Node& Tree::GetNode(size_t i) const {
	return linearizedNodes_[i];
}

void Tree::print(ostream& os) const {
	for (auto it = this->rbegin(); it !=  this->rend(); it++) {
		const Node& node = *it;
		node.info(os);
		node.shape().print(os);
		os << endl;
	}
}

ostream& operator<<(ostream& os, const Tree& basis) {
	if(&os == &cout) {
		basis.info(os);
	} else {
		basis.Write(os);
	}
	return os;
}

istream& operator>>(istream& is, Tree& basis) {
	basis.Read(is);
	return is;
}

