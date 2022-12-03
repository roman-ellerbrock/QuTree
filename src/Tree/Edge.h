//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef EDGE_H
#define EDGE_H
#include "Tree/Node.h"


class Edge {
	/**
	 * \class Edge
	 * \brief This class describes the edge between two nodes.
	 * 
	 * Example 1:
	 * Node x [from] (child, Childidx: 2, parentIdx: 3) ---> Node y [to] (parent)
	 *
	 * inIdx: incoming to *to_* over index (parentIdx) 3
	 * outIdx: outgoing to *to_* over index (childIdx) 2
	 *
	 * Example 2:
	 * Node x [to] (child, Childidx: 2, parentIdx: 3) <--- Node y [from] (parent)
	 *
	 * inIdx: incoming to *to_* over index (childIdx) 2
	 * outIdx: outgoing to *to_* over index (parentIdx) 3
	 */
public:
	Edge() = default;

	Edge(const Node& from, const Node& to)
		: from_(&from), to_(&to) {}

	Edge(const Node* from, const Node* to)
		: from_(from), to_(to) {}

	[[nodiscard]] const Node& from() const {
		if (!from_) {
			cerr << "Edge error: cannot convert nullptr (from_) to Node&.\n";
			exit(1);
		}
		return *from_;
	}

	[[nodiscard]] const Node& to() const {
		if (!to_) {
			cerr << "Edge error: cannot convert nullptr (to_) to Node&.\n";
			exit(1);
		}
		return *to_;
	}

	/// index over which *from_* goes out to *to_*
	[[nodiscard]] size_t outIdx() const {
		if (!from_) {
			cerr << "Edge error: cannot evaluate outIdx (from_ is nullptr).\n";
			exit(1);
		}
		return from_->idx(*to_);
	}

	/// index over which *from_* comes in to *to_*
	[[nodiscard]] size_t inIdx() const {
		if (!from_) {
			cerr << "Edge error: cannot evaluate inIdx (to_ is nullptr).\n";
			exit(1);
		}
		return to_->idx(*from_);
	}

	[[nodiscard]] bool isUpEdge() const {
		if (from_->parent_ == to_) {
			return true;
		} else {
			return false;
		}
	}

	[[nodiscard]] const Node& up() const {
		if (isUpEdge()) {
			return to();
		} else  {
			return from();
		}
	}

	[[nodiscard]] const Node& down() const {
		if (isUpEdge()) {
			return from();
		} else  {
			return to();
		}
	}

	bool operator==(const Edge& b) const {
		if (from() != b.from()) { return false; }
		if (to() != b.to()) { return false; }
		return true;
	}

	bool operator!=(const Edge& b) const {
		return !(*this == b);
	}

	/*
	 * Note: the address of edges does not reflect it's position in a sweep or
	 *       in a vector of edge-attributes! We will change this in the future.
	 *       The address is built in a way that it doesn't require of the knowledge
	 *       of how many nodes there are in the graph.
	 */
	[[nodiscard]] size_t address() const {
		return 2 * down().address_ - (!isUpEdge()) - 1;
	}

	// This is a uid for the vector of down edges OR up_edges
	[[nodiscard]] size_t local_address() const {
		return down().address_ - 1;
	}

	[[nodiscard]] TensorShape shape() const {
		size_t dim = from().shape_[outIdx()];
		return TensorShape({dim, dim});
	}

	void info(ostream& os = cout) const {
		from_->position_.info(os, false, string(" "));
		os << "--> ";
		to_->position_.info(os, false, string(" "));
		os << "(id: " << address() << ")\n";
	}

private:
	const Node *from_;
	const Node *to_;
};

Edge outgoingEdge(const Node& from, size_t toIdx);
Edge incomingEdge(const Node& to, size_t inIdx);

vector<Edge> outgoingEdges(const Node& from);
vector<Edge> incomingEdges(const Node& to);

vector<Edge> preEdges(const Edge& edge);
vector<Edge> preUpEdges(const Edge& e);
vector<Edge> preDownEdges(const Edge& e);

vector<Edge> postEdges(const Edge& edge);
vector<Edge> postEdges(const Edge* edge);

ostream& operator<<(ostream& os, const Edge& edge);

bool operator==(const vector<Edge>& a, const vector<Edge>& b);
bool operator!=(const vector<Edge>& a, const vector<Edge>& b);

ostream& operator<<(ostream& os, const vector<Edge>& edges);

#endif //EDGE_H
