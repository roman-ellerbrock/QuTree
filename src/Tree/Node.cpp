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
	  address_(node.address_) {
	reconnect();
}

Node& Node::operator=(const Node& old) {
	Node node(old);
	swap(*this, node);
	reconnect();
	return *this;
}

Node& Node::operator=(Node&& old) noexcept {
	if (this == &old) {
		return *this;
	}
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

	read(file, up, position_);
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

Node read(istream& file, Node *up,
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
			node.child_.emplace_back(read(file, &node, newposition));
		} else {
			node.leaves_.emplace_back(Leaf(file, &node, newposition));
		}
	}

	// create a TensorDim after dimensions were read
	dim.push_back(nstates);
	node.shape_ = TensorShape(dim);
	return node;
}

/*Node::Node(const Leaf& leaf, size_t ntensor)
	: Node() {
	leaves_.push_back(leaf);
	reconnect();

	shape_ = TensorShape({leaf.api_.basis()->par_.dim_, ntensor});
}*/

/*Node *Node::nextNode() {

	Node *result;
	if (nextNodeNum_ >= 0) {
		result = child_[nextNodeNum_].nextNode();
		if (result == &child_[nextNodeNum_]) {
			nextNodeNum_--;
		}
	} else {
		nextNodeNum_ = child_.size() - 1;
		result = this;
	}

	return result;
}

Node *Node::nextNodeManthe() {

	Node *result;
	if (nextNodeNumFortran_ < child_.size()) {
		result = child_[nextNodeNumFortran_].nextNodeManthe();
		if (result == &child_[nextNodeNumFortran_]) {
			++nextNodeNumFortran_;
		}
	} else {
		nextNodeNumFortran_ = 0;
		result = this;
	}

	return result;
}*/


/*void Node::expandChild(size_t i) {
	assert(!isBottomlayer());
	assert(i < child_.size());

	// Make a new "down_new" by adding all old children
	// except the expanded node and adding the expanded node's children
	// at the right place

	vector<unique_ptr<AbstractNode>> down_new;
	// move nodes until expanded node to down_new
	for (size_t j = 0; j < i; j++) {
		down_new.push_back(move(child_[j]));
	}

	// append the children of the expanded node to down_new
	Node& child = this->child(i);
	assert(!child.isBottomlayer());
	size_t nchildren = child.nChildren();
	for (size_t j = 0; j < nchildren; j++) {
		// Update parent-ptr of new children
		Node& subnode = child.child(j);
		subnode.setParent(this);
		// Update subnodes positionindices
		int child_nr = down_new.size();
		subnode.updatePosition(position_ * child_nr);

		// Save the unique_ptr to the new vector
		unique_ptr<AbstractNode> subchild = child.downUnique(j);
		down_new.push_back(move(subchild));
	}

	// move the rest of the old children to down_new
	for (size_t j = i + 1; j < nChildren(); j++) {
		down_new.push_back(move(child_[j]));
	}

	// move down_new to down_
	child_ = move(down_new);

	// Adjust TensorDim
	updateTDim();

	// Update position
	updatePosition(position_);

	// update nTotalNodes_: All nTotalNodes_ of the nodes above change
	// so go to the topnode and update from there downwards
	Node& topnode = topNode();
	topnode.updatennodes();

	// Set the nextNodeNum_ for a correct nextNode() sweep that starts
	// at the last_ node
	nextNodeNum_ = child_.size() - 1;
}
*/

