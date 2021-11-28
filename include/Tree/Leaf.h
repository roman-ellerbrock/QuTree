#pragma once
#include "stdafx.h"
#include "NodePosition.h"
#include "PrimitiveBasis/LeafAPI.h"

class Node;

class Leaf
	/**
	 * \class Leaf
	 * \ingroup Tree
	 * \brief This class represents the dangling edges in a tree.
	 *
	 * A Leaf is a dangling edge of a TTN in a Tree.
	 * in a tree of a TTBasis. Leaves contain abstract class pointers
	 * to PrimitiveBasis which provides the interface to the problem
	 * under consideration.
	 */
{
public:
	Leaf()
		: parent_(nullptr), position_(), api_() {};
	~Leaf() = default;

	Leaf(istream& is, Node* parent, const NodePosition& pos) {
		//@TODO: add test!
		parent_ = parent;
		position_ = pos;

		BasisParameters par;
		par.readDimline(is);
		initialize(par);
	}

	void initialize(const BasisParameters& par) {
		api_.initialize(par);
	}

	void info(ostream& os = cout) const { position_.info(os); }

	void write(ostream& os) const {
		for (size_t l = 0; l < position_.layer(); l++) { os << "\t"; }
		api_.write(os);
	}

	BasisParameters& par() {
		assert(api_.basis());
		return api_.basis()->par_;
	}

	[[nodiscard]] const BasisParameters& par() const {
		assert(api_.basis());
		return api_.basis()->par_;
	}

	[[nodiscard]] const Node& parent() const { return *parent_; };

	bool operator==(const Leaf& b) const {
		if (this->position_ != b.position_) { return false; }
		if (this->par() != b.par()) { return false; }
		return true;
	}

	bool operator!=(const Leaf& b) const {
		return !(*this == b);
	}

	Node *parent_;
	NodePosition position_;
	LeafAPI api_;
};

