//
// Created by Roman Ellerbrock on 11/30/21.
//

#include <Tree/Edge.h>

Edge outgoingEdge(const Node& from, size_t outIdx) {
	if (outIdx < from.nChildren()) {
		const Node& to = from.child_[outIdx];
		return {from, to};
	} else if (outIdx == from.parentIdx()){
		const Node& to = from.parent();
		return {from, to};
	} else {
		cerr << "Error creating outgoing edge.\n";
		exit(1);
	}
}

Edge incomingEdge(const Node& to, size_t inIdx) {
	if (inIdx < to.nChildren()) {
		const Node& from = to.child_[inIdx];
		return {from, to};
	} else if (inIdx == to.parentIdx()){
		const Node& from = to.parent();
		return {from, to};
	} else {
		cerr << "Error creating incoming edge.\n";
		exit(1);
	}
}

vector<Edge> outgoingEdges(const Node& from) {
	vector<Edge> edges;
	for (const Node& child : from.child_) {
		edges.emplace_back(Edge(from, child));
	}
	if (from.parent_) {
		edges.emplace_back(Edge(from, from.parent()));
	}
	return edges;
}

vector<Edge> incomingEdges(const Node& to) {
	vector<Edge> edges;
	for (const Node& child : to.child_) {
		edges.emplace_back(Edge(child, to));
	}
	if (to.parent_) {
		edges.emplace_back(Edge(to.parent(), to));
	}
	return edges;
}

vector<Edge> preUpEdges(const Edge& e) {
	const Node& to = e.from();
	vector<Edge> edges;
	for (const Node& from : to.child_) {
		edges.emplace_back(Edge(from, to));
	}
	return edges;
}

vector<Edge> preDownEdges(const Edge& e) {
	const Node& center = e.from();
	size_t hole = e.outIdx();

	vector<Edge> edges;
	for (size_t k = 0; k < center.nChildren(); ++k) {
		if (k == hole) { continue; }
		const Node& from = center.child_[k];
		edges.emplace_back(Edge(from, center));
	}
	if (center.parent_) {
		const Node& parent = center.parent();
		edges.emplace_back(Edge(parent, center));
	}
	return edges;
}

vector<Edge> preEdges(const Edge& edge) {
	if (edge.isUpEdge()) {
		return preUpEdges(edge);
	} else {
		return preDownEdges(edge);
	}
}

ostream& operator<<(ostream& os, const Edge& edge) {
	os << "Edge(" << edge.from().address_ <<", " << edge.to().address_ << ")";
	return os;
}

bool operator==(const vector<Edge>& a, const vector<Edge>& b) {
	if (a.size() != b.size()) { return false; }
	for (size_t i = 0; i < a.size(); ++i) {
		if (a[i] != b[i]) { return false; }
	}
	return true;
}

ostream& operator<<(ostream& os, const vector<Edge>& edges) {
	for (const Edge& e : edges) {
		cout << e << " ";
	}
	cout << endl;
	return os;
}

