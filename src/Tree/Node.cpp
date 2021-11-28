#include "Tree/Node.h"

Node::Node()
	: shape_({1}),
	  address_(-100),
	  parent_(nullptr),
	  position_() {}

Node::Node(const Node& node)
	: parent_(node.parent_),
	  position_(node.position_),
	  address_(node.address_),
	  child_(node.child_),
	  leaves_(node.leaves_),
	  shape_(node.shape_) {

	reconnect();
}

Node::Node(Node&& node) noexcept
	: parent_(node.parent_),
	  child_(move(node.child_)),
	  leaves_(move(node.leaves_)),
	  position_(node.position_),
	  shape_(node.shape_),
	  address_(node.address_) {
	reconnect();
}

Node& Node::operator=(const Node& other) {
	if (this == &other) {
		return *this;
	}

	parent_ = other.parent_;
	child_ = other.child_;
	leaves_ = other.leaves_;
	position_ = other.position_;
	shape_ = other.shape_;
	address_ = other.address_;
	reconnect();

	return *this;
}

Node& Node::operator=(Node&& old) noexcept {
	swap(*this, old);
	reconnect();
	return *this;
}

void Node::reconnect() {
	size_t k = 0;
	for (Node& child : child_) {
		child.parent_ = this;
		child.updatePosition(position_ * k++);
	}
	for (Leaf& leaf : leaves_) {
		leaf.parent_ = this;
		leaf.position_ = position_ * k++;
	}
}

Node::Node(istream& file, Node *up,
	const NodePosition& position)
	: parent_(up),
	  position_(position),
	  address_(-100) {

	readNode(file, up, position_);
}

void Node::info(ostream& os) const {
	position_.info(os);
//	shape_.print(os);
}

void Node::write(ostream& os) const {
	const TensorShape& tdim = shape_;
	for (size_t l = 0; l < position_.layer(); l++) { os << "\t"; }
	os << tdim.lastDimension() << "\t-" << nChildren() << "\n";
	for (const Node& child: child_) {
		child.write(os);
	}
}

Node& Node::parent() {
	assert(parent_ != nullptr);
	return *parent_;
}

const Node& Node::parent() const {
	assert(parent_ != nullptr);
	return *parent_;
}

const Node *Node::get(NodePosition p) const {
	if (p.empty()) {
		return this;
	} else {
		size_t k = p.front();
		p.pop_front();
		assert(k < child_.size());
		return child_[k].get(p);
	}
}

void Node::updatePosition(const NodePosition& p) {
	// Save the new position_
	position_ = p;

	size_t k = 0;
	for (; k < nChildren(); ++k) {
		NodePosition subp = p * k;
		child_[k].updatePosition(subp);
	}
	for (Leaf& leaf : leaves_) {
		NodePosition subp = p * k++;
		leaf.position_ = subp;
	}
}

Node *askParent(const Node *p) {
	if (p->parent_ == nullptr) { return nullptr; }
	size_t next = p->childIdx() + 1;
	if (next < p->parent_->nChildren()) {
		return &p->parent_->child_[next];
	} else {
		return askParent(p->parent_);
	}
}

Node *sweep(Node *p) {
	if (!p->child_.empty()) {
		return &(p->child_.front());
	} else {
		return askParent(p);
	}
}

Node readNode(istream& file, Node *up,
	const NodePosition& position) {
	Node node;
	node.parent_ = up;
	node.position_ = position;

	// Variables needed to read input file
	int nstates, f;
	vector<size_t> dim;

	// read input file
	file >> nstates;
	assert(nstates > 0);
	file >> f;
	f = -f;
	assert(f > 0);

	for (int i = 0; i < f; i++) {
		// check wether physical or logical coord
		// and reset stream to old position_
		int newstate, newf;
		auto mark = file.tellg();
		file >> newstate;
		assert(newstate > 0);
		file >> newf;
		file.seekg(mark);

		dim.push_back(newstate);

		// generate position_ for the new node
		NodePosition newposition = node.position_ * i;

		if (newf < 0) {
			node.child_.emplace_back(readNode(file, &node, newposition));
		} else {
			node.leaves_.emplace_back(Leaf(file, &node, newposition));
		}
	}

	// create a TensorDim after dimensions were read
	dim.push_back(nstates);
	node.shape_ = TensorShape(dim);
	return node;
}

bool operator==(const Node& a, const Node& b) {
	/// only test for parent equality if a||b is nullptr
	if ((a.parent_== nullptr) && (b.parent_!=nullptr)) { return false; }
	if ((b.parent_== nullptr) && (a.parent_!=nullptr)) { return false; }
	if (a.shape_ != b.shape_) { return false; }
	if (a.address_ != b.address_) { return false; }
	if (a.position_ != b.position_) { return false; }
	if (a.child_ != b.child_) { return false; }
	if (a.leaves_ != b.leaves_) { return false; }
	return true;
}

bool operator!=(const Node& a, const Node& b) {
	return !(a == b);
}
