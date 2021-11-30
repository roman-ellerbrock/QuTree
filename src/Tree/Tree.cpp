//
// Created by Roman Ellerbrock on 2020-01-21.
//
#include "Tree/Tree.h"

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

void Tree::update() {
	root_.updatePosition(NodePosition());
	nodeArray_ = NodeArray(root_);
	leafArray_ = LeafArray(root_);
	edgeArray_ = EdgeArray(nodeArray_);
}

void Tree::read(istream& is) {
	root_ = readNode(is);
	update();
	leafArray_.readPars(is);
}

void Tree::read(const string& filename) {
	ifstream is(filename);
	read(is);
}

void Tree::write(ostream& os) const {
	root_.write(os);
	leafArray_.writePars(os);
}

void Tree::info(ostream& os) const {
	os << "List of Leaves:" << endl;
	for (size_t i = 0; i < this->nLeaves(); i++) {
		const Leaf& leaf = leafArray_[i];
		leaf.info(os);
		os << endl;
	}
	os << endl;

	os << "List of upper nodes:" << endl;
	for (int i = nNodes() - 1; i >= 0; i--) {
		const Node& node = nodeArray_[i];
		node.info();
		node.shape_.print(os);
		os << endl;
	}
	os << "Number of Nodes = " << nNodes() << endl;
}

void Tree::print(ostream& os) const {
	for (auto it = this->rbegin(); it != this->rend(); it++) {
		const Node& node = *it;
		node.info(os);
		node.shape_.print(os);
		os << endl;
	}
}

ostream& operator<<(ostream& os, Tree& tree) {
	tree.write(os);
	return os;
}

istream& operator<<(istream& is, Tree& tree) {
	tree.read(is);
	return is;
}

ostream& operator<<(ostream& os, const Tree& tree) {
	if (&os == &cout) {
		tree.info(os);
	} else {
		tree.write(os);
	}
	return os;
}

istream& operator>>(istream& is, Tree& tree) {
	tree.read(is);
	return is;
}

bool operator==(const Tree& a, const Tree& b) {
	if (a.root() != b.root()) { return false; }
	if (a.leafArray().size() != b.leafArray().size()) { return false; }
	if (a.nodeArray().size() != b.nodeArray().size()) { return false; }
	return true;
}



/*void Tree::replaceNode(Node& old_node, Node& new_node) {
	// The old node must not be the toplayer node, otherwise change the
	// whole tree
	assert(!old_node.isToplayer());

	// Replace the node
	Node& parent = old_node.parent();
	parent.replace(new_node, old_node.childIdx());

	Node& topnode = topNode();
	topnode.updatePosition(NodePosition());
	linearizeLeaves();
	linearizeNodes();
}*/

/*bool Tree::isWorking() {

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
}*/
/*void Tree::resetLeafModes() {
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
}*/

/*void Tree::reindexLeafModes(map<size_t, size_t> Map) {
	for (Node& node : *this) {
		if (node.isBottomlayer()) {
			Leaf& leaf = node.getLeaf();
			leaf.mode() = Map[leaf.mode()];
		}
	}
	update();
}*/

