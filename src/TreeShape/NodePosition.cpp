#include "TreeShape/NodePosition.h"

void NodePosition::info(ostream& os, bool print_layer) const {
	if (print_layer) {
	os << "{ " << layer_ << " ; ";
	} else {
		os << "{ ";
	}
	for (size_t i = 0; i < path_.size(); i++) {
		if (i > 0) { os << ", "; }
		os << path_[i];
	}
	os << " }" << endl;
}

size_t NodePosition::childIdx() const {
	size_t size = path_.size();
	if (size > 0) {
		return path_[size - 1];
	} else {
		return 0;
	}
}

NodePosition operator*(NodePosition p, size_t k) {
	// contruct p_+k by going from p_ to direction k
	NodePosition pnew(p);
	pnew.layer_++;
	pnew.push_back(k);
	return pnew;
}

NodePosition operator*(NodePosition p, NodePosition q) {
	// contruct p_+q by going from root to p_ and then to q
	NodePosition pnew(p);
	pnew.layer_ = p.layer_ + q.layer_;
	for (size_t i = 0; i < q.path_.size(); i++) {
		pnew.push_back(q.path_[i]);
	}
	return pnew;
}

