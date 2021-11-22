#include "TreeShape/Node.h"

Node::Node()
	: nodeType_(1), bottomLayer_(false), nextNodeNum_(0),
	nextNodeNumFortran_(0), address_(-100), parent_(nullptr) {
}

// Copy constructor
Node::Node(const Node& node)
	: nTotalNodes_(node.nTotalNodes_), nNodes_(node.nNodes_),
	  nLeaves_(node.nLeaves_), tensorDim_(node.tensorDim_),
	  parent_(node.parent_),
	  nextNodeNum_(node.nextNodeNum_),
	  nextNodeNumFortran_(node.nextNodeNumFortran_),
	  position_(node.position_),
	  address_(node.address_), nodeType_(node.nodeType_),
	  bottomLayer_(node.bottomLayer_) {
	if (isBottomlayer()) {
		child_.emplace_back(make_unique<Leaf>(node.getLeaf()));
		getLeaf().setParent(this);
	} else {
		for (size_t i = 0; i < node.nChildren(); i++) {
			child_.emplace_back(make_unique<Node>(node.child(i)));
			child(i).setParent(this);
		}
	}
}

// move constructor
Node::Node(Node&& node) noexcept
	: nTotalNodes_(node.nTotalNodes_), nNodes_(node.nNodes_),
	  nLeaves_(node.nLeaves_), tensorDim_(node.tensorDim_),
	  parent_(node.parent_),
	  nextNodeNum_(node.nextNodeNum_), position_(node.position_),
	  address_(node.address_), nodeType_(node.nodeType_),
	  bottomLayer_(node.bottomLayer_), nextNodeNumFortran_(node.nextNodeNumFortran_) {
	// Set upwards connectivity of children
	child_ = move(node.child_);
	if (node.isBottomlayer()) {
		getLeaf().setParent(this);
	} else {
		for (size_t i = 0; i < nChildren(); i++) {
			child(i).setParent(this);
		}
	}
}

Node& Node::operator=(const Node& old) {
	Node node(old);
	*this = move(node);
	return *this;
}

Node& Node::operator=(Node&& old) noexcept {
	if (this == &old) {
		return *this;
	}
	nTotalNodes_ = old.nTotalNodes_;
	nNodes_ = old.nNodes_;
	nLeaves_ = old.nLeaves_;
	tensorDim_ = old.tensorDim_;
	parent_ = old.parent_;
	nextNodeNum_ = old.nextNodeNum_;
	position_ = old.position_;
	address_ = old.address_;
	nodeType_ = old.nodeType_;
	bottomLayer_ = old.bottomLayer_;
	child_ = move(old.child_);
	// Set upwards connectivity of children
	if (isBottomlayer()) {
		getLeaf().setParent(this);
	} else {
		for (size_t i = 0; i < nChildren(); i++) {
			child(i).setParent(this);
		}
	}

	return *this;
}

Node::Node(const Leaf& phys, size_t ntensor)
	: Node() {
	// This constructor initializes a bottomLayer_ GetNode from a
	// Leaf
//	down_.emplace_back(unique_ptr<Leaf>(new Leaf(phys)));
	child_.emplace_back(make_unique<Leaf>(phys));
	nTotalNodes_++;
	nLeaves_++;
	bottomLayer_ = true;
	nextNodeNum_ = 0;
	nextNodeNumFortran_ = 0;

	// Set connectivity of linearizedLeaves_ node
	Leaf& phy = getLeaf();
	phy.setParent(this);

	// Build the TensorDim
	tensorDim_ = TensorShape({phys.dim(), ntensor});
}

void Node::initialize(istream& file, Node *up,
	const NodePosition& position) {
	parent_ = up;
	position_ = position;

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
		NodePosition newposition = position_ * i;

		if (newf < 0) {
			// add logical node
			child_.emplace_back(unique_ptr<Node>
				(new Node(file, this, newposition)));
			nTotalNodes_ += (*child_[i]).nTotalNodes();
			nNodes_ += (*child_[i]).nNodes();
			nLeaves_ += (*child_[i]).nLeaves();
		} else {
			// add physical node
			child_.emplace_back(unique_ptr<Leaf>
				(new Leaf(file, this, newposition)));
			// A physical node was added, so increase the number of nodes by 1
			nTotalNodes_++;
			nLeaves_++;
			// If this layer_ holds a Physical Coordinate it is a bottom layer_
			bottomLayer_ = true;
		}
	}

	// create a TensorDim after dimensions were read
	dim.push_back(nstates);
	tensorDim_ = TensorShape(dim);

	nextNodeNum_ = child_.size() - 1;
}

Node::Node(istream& file, Node *up,
	const NodePosition& position)
	: parent_(up), nTotalNodes_(1),
	  nLeaves_(0), nNodes_(1), position_(position),
	  nodeType_(1), bottomLayer_(false), nextNodeNum_(0),
	  nextNodeNumFortran_(0), address_(-100) {
	// Call Initialize
	initialize(file, up, position_);
}

void Node::info(ostream& os) const {
	position_.info(os);
//	tensorDim_.print(os);
}

void Node::write(ostream& os) const {
	const TensorShape& tdim = shape();
	for (size_t l = 0; l < position_.layer(); l++) { os << "\t"; }
	os << tdim.lastDimension() << "\t-" << nChildren() << "\n";
	for (size_t i = 0; i < nChildren(); i++) {
		child_[i]->write(os);
	}
}

void Node::push_back(const Node& node) {
	child_.emplace_back(std::make_unique<Node>(node));
	child_.back()->setParent(this);
	updatennodes();
}

bool Node::isToplayer() const {
	return (parent_ == nullptr);
}

Leaf& Node::getLeaf() {
	assert(bottomLayer_);
	AbstractNode *node = child_[0].get();
	return (Leaf&) (*node);
}

const Leaf& Node::getLeaf() const {
	assert(bottomLayer_);
	AbstractNode *node = child_[0].get();
	return (Leaf&) (*node);
}

const Node& Node::child(size_t i) const {
	assert(i < child_.size());
	assert(!bottomLayer_);
	return (Node&) *child_[i];
}

Node& Node::child(size_t i) {
	assert(i < child_.size());
	assert(!bottomLayer_);
	return (Node&) *child_[i];
}

Node& Node::parent() {
	assert(parent_ != nullptr);
	return (Node&) *parent_;
}

const Node& Node::parent() const {
	assert(parent_ != nullptr);
	return (Node&) *parent_;
}

AbstractNode *Node::nextNode() {

	AbstractNode *result;
	if (nextNodeNum_ >= 0) {
		result = child_[nextNodeNum_].get()->nextNode();
		if (result == child_[nextNodeNum_].get()) {
			nextNodeNum_--;
		}
	} else {
		nextNodeNum_ = child_.size() - 1;
		result = this;
	}

	return result;
}

AbstractNode *Node::nextNodeManthe() {

	AbstractNode *result;
	if (nextNodeNumFortran_ < child_.size()) {
		result = child_[nextNodeNumFortran_].get()->nextNodeManthe();
		if (result == child_[nextNodeNumFortran_].get()) {
			++nextNodeNumFortran_;
		}
	} else {
		nextNodeNumFortran_ = 0;
		result = this;
	}

	return result;
}

unique_ptr<AbstractNode> Node::downUnique(size_t i) {
	assert(i < child_.size());
	return move(child_[i]);
}

void Node::expandChild(size_t i) {
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

void Node::update(const NodePosition& p) {
	// @TODO: Should reset state_ and Update(connectivity) be in separate routines?
	resetCounters();
	updatePosition(p);
	updatennodes();
	updateTDim();
}

void Node::updateTDim() {
/*	size_t ntensor = 0;
	if (IsToplayer()) {
		// If the node is a Toplayer, ntensor stays the same
		ntensor = tensorDim_.getntensor();
	} else {
		// If this node is not a toplayer-node, ntensor is given by the parents 
		// active_ size
		Node& parent = parent();
		// @TODO: This looks wrong - check again. Doesnt it have to be active_(k)?
		ntensor = parent.shape().LastActive();
	}*/

	// Get the dimensions of the children by requesting their ntensors
	vector<size_t> dim_new;
	if (isBottomlayer()) {
		const Leaf& phys = getLeaf();
		dim_new.push_back(phys.dim());
	} else {
		for (int i = 0; i < nChildren(); i++) {
			const Node& child = this->child(i);
			const TensorShape& tdimchild = child.shape();
			dim_new.push_back(tdimchild.lastDimension());
		}
	}

	// Create a new TensorDim from the dim_-vector and ntensor
	dim_new.push_back(tensorDim_.lastDimension());
	tensorDim_ = TensorShape(dim_new);
}

void Node::updatePosition(const NodePosition& p) {
	// Save the new position_
	position_ = p;

	// Update positions of children
	for (size_t i = 0; i < child_.size(); i++) {
		NodePosition subp = p * i;
		if (isBottomlayer()) {
			Leaf& phy = getLeaf();
			phy.updatePosition(subp);
		} else {
			Node& child = this->child(i);
			child.updatePosition(subp);
		}
	}
}

void Node::resetCounters() {

	nextNodeNum_ = child_.size() - 1;
	nextNodeNumFortran_ = 0;
	if (!isBottomlayer()) {
		for (size_t k = 0; k < nChildren(); ++k) {
			child(k).resetCounters();
		}
	}
}

void Node::updatennodes() {
	if (isBottomlayer()) {
		nTotalNodes_ = 2;
		nNodes_ = 1;
		nLeaves_ = 1;
	} else {
		nTotalNodes_ = 1;
		nNodes_ = 1;
		nLeaves_ = 0;
		for (size_t k = 0; k < child_.size(); k++) {
			Node& child = this->child(k);
			child.updatennodes();
			nTotalNodes_ += child.nTotalNodes();
			nNodes_ += child.nNodes();
			nLeaves_ += child.nLeaves();
		}
	}
}

Node& Node::topNode() {
	// Returns the topnode of the tree
	if (isToplayer()) {
		return (*this);
	} else {
		return parent().topNode();
	}
}

void Node::replace(Node& new_child, size_t idx) {
	assert(idx < child_.size());
	child_[idx] = unique_ptr<Node>(new Node(new_child));
}



