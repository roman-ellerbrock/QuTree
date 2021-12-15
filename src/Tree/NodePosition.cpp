#include "Tree/NodePosition.h"

void NodePosition::info(ostream& os, bool print_layer, const string& term) const {
	if (print_layer) {
	os << "{ " << layer() << " ; ";
	} else {
		if (empty()) {
			os << "{ r";
		} else {
			os << "{ r, ";
		}
	}
	size_t i = 0;
	for (auto x : *this) {
		if (i++ > 0) { os << ", "; }
		os << x;
	}
	os << " }" << term;
}

size_t NodePosition::childIdx() const {
	assert(!empty());
	return back();
}

NodePosition operator*(NodePosition p, size_t k) {
	// contruct p_+k by going from p_ to direction k
	NodePosition pnew(p);
	pnew.push_back(k);
	return pnew;
}

NodePosition operator*(NodePosition p, NodePosition q) {
	// contruct p_+q by going from root to p_ and then to q
	p.splice(p.end(), q);
	return p;
}

bool operator==(NodePosition p, NodePosition q) {
	if (p.size() != q.size()) { return false; }
	for (size_t i = 0; i < p.size(); ++i) {
		size_t x = p.back();
		size_t y = q.back();
		p.pop_back();
		q.pop_back();
		if (x != y) { return false; }
	}
	return true;
}

ostream& operator<<(ostream& stream, const NodePosition& p) {
	p.info(stream);
	return stream;
}
