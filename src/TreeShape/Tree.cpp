//
// Created by Roman Ellerbrock on 2020-01-21.
//
#include "TreeShape/Tree.h"

Tree::Tree(const Tree& T)
	: root_(T.root_) {
	update();
}

Tree::Tree(Tree&& T) noexcept {
	root_ = move(T.root_);
	update();
}

Tree& Tree::operator=(const Tree& T) {
	*this = Tree(T);
	return *this;
}

Tree& Tree::operator=(Tree&& T) noexcept {
	root_ = move(T.root_);
	update();
	return *this;
}

Tree::Tree(const string& filename) {
	read(filename);
}

Tree::Tree(istream& is) {
	read(is);
}

void Tree::resetLeafModes() {
	size_t n_modes = this->nLeaves();
	assert(n_modes > 0);
	int mode = n_modes - 1;
	for (Node& node : *this) {
		if (node.isBottomlayer()) {
			Leaf& leaf = node.getLeaf();
			leaf.mode() = mode--;
		}
	}
	update();
}

void Tree::reindexLeafModes(map<size_t, size_t> Map) {
	for (Node& node : *this) {
		if (node.isBottomlayer()) {
			Leaf& leaf = node.getLeaf();
			leaf.mode() = Map[leaf.mode()];
		}
	}
	update();
}

void Tree::expandNode(Node& node) {
	assert(!node.isToplayer());
	assert(!node.isBottomlayer());

	Node& parent = node.parent();
	size_t childIdx = node.childIdx();
	parent.expandChild(childIdx);
	linearizeNodes();
}

void Tree::update() {
	// Tree is assumed to be updated, but the rest not:
	// Update everything
	root_.update(NodePosition());
	linearizeLeaves();
	linearizeNodes();
}

void Tree::replaceNode(Node& old_node, Node& new_node) {
	// The old node must not be the toplayer node, otherwise change the
	// whole tree
	assert(!old_node.isToplayer());

	// Replace the node
	Node& parent = old_node.parent();
	parent.replace(new_node, old_node.childIdx());

	Node& topnode = topNode();
	topnode.updatePosition(NodePosition());
	topnode.updatennodes();
	linearizeLeaves();
	linearizeNodes();
}

void Tree::linearizeNodes() {
	// block has to be cleared, because logical block
	// must be resistant to re-feed (important for e.g. expand node)
	// This routine adds every mctdh node to the logical block
	linearizedNodes_.clear();
	int counter = 0;
	for (int i = 0; i < nTotalNodes(); i++) {
		AbstractNode& abstract_node = nextNode();
		if (abstract_node.type() == 1) {
			auto& node = (Node&) (abstract_node);
			node.setAddress(counter);
			counter++;
			linearizedNodes_.push_back(node);
		}
	}

	edges_.clear();
	for (const Node& node : linearizedNodes_) {
		if (!node.isToplayer()) {
			const Node& parent = node.parent();
			edges_.emplace_back(Edge(node, parent));
		}
	}
}

void Tree::linearizeLeaves() {
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

void Tree::read(istream& is) {
	// feed linearizedLeaves_ and logical block with references
	root_.initialize(is, nullptr, NodePosition());
	update();

	// Add new PhysPar for every physical coordinate
	for (int i = 0; i < linearizedLeaves_.size(); i++) {
		// Set parameters
		PhysPar par(is);
		linearizedLeaves_[i].setPar(par);

		// Initialize primitive grid (HO, FFT, Legendre, ...)
		LeafInterface& primitivebasis = linearizedLeaves_[i].interface();
		primitivebasis.initialize(par.omega(), par.r0(), par.wfr0(), par.wfOmega());
	}
}

void Tree::read(const string& filename) {
	ifstream is(filename);
	read(is);
}

void Tree::write(ostream& os) const {
	root_.write(os);
}

ostream& operator<<(ostream& os, Tree& basis) {
	basis.write(os);
	return os;
}

istream& operator<<(istream& is, Tree& basis) {
	basis.read(is);
	return is;
}

void Tree::info(ostream& os) const {
	os << "List of Leaves:" << endl;
	for (size_t i = 0; i < this->nLeaves(); i++) {
		const Leaf& node = getLeaf(i);
		node.info(os);
		os << endl;
	}
	os << endl;

	// ... and now for every logical node
	os << "List of upper nodes:" << endl;
	for (int i = nNodes() - 1; i >= 0; i--){
		const Node& node = getNode(i);
		node.info();
		node.shape().print(os);
		os << endl;
	}
	os << "Number of Nodes = " << nNodes() << endl;
}

bool Tree::isWorking() {

	bool works = true;
	int counter = 0;
	for (int i = 0; i < nTotalNodes(); i++) {
		const AbstractNode& abstract_node = nextNode();
		if (abstract_node.type() == 1) {
			auto& node = (Node&) (abstract_node);
			// 1.) Check global address
			// Do not break here to leave nodes in a valid state
			if (counter != node.address()) { works = false; }
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
		cerr << "Corrupted linearizedNodes_. Missing update()?" << endl;
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
			if (&leaf != &linearizedLeaves_[leaf.mode()]) { works = false; }
		}
	}
	if (!works) {
		cerr << "Corrupted linearizedLeaves." << endl;
		return false;
	}

	return true;
}

Leaf& Tree::getLeaf(size_t i) {
	return linearizedLeaves_[i];
}

const Leaf& Tree::getLeaf(size_t i) const {
	return linearizedLeaves_[i];
}

Node& Tree::getNode(size_t i) {
	return linearizedNodes_[i];
}

const Node& Tree::getNode(size_t i) const {
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
		basis.write(os);
	}
	return os;
}

istream& operator>>(istream& is, Tree& basis) {
	basis.read(is);
	return is;
}

vector<const Node*> Tree::neighbors(const Node& from, int hole) const {
	vector<const Node*> neigh;
	for (size_t k = 0; k < from.nChildren(); ++k) {
		if (k == hole) { continue; }
		neigh.push_back(&from.child(k));
	}

	if (!from.isToplayer()) {
		if (from.parentIdx() != hole) {
			neigh.push_back(&from.parent());
		}
	}
	return neigh;
}
