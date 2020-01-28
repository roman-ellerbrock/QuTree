#include "TensorTreeBasis/NodePosition.h"

void NodePosition::info(ostream& os) const {
	os << "{ " << layer_ << " ; ";
	for (int i = 0; i < path_.size(); i++) {
		if (i > 0) { os << ", "; }
		os << path_[i];
	}
	os << " }" << endl;
}

int NodePosition::ChildIdx() const {
	int size = path_.size();
	if (size > 0) {
		return path_[size - 1];
	} else {
		return 0;
	}
}

NodePosition operator*(NodePosition p, int k) {
	// contruct p+k by going from p to direction k
	NodePosition pnew(p);
	pnew.layer_++;
	pnew.push_back(k);
	return pnew;
}

NodePosition operator*(NodePosition p, NodePosition q) {
	// contruct p+q by going from root to p and then to q
	NodePosition pnew(p);
	pnew.layer_ = p.layer_ + q.layer_;
	for (int i = 0; i < q.path_.size(); i++) {
		pnew.push_back(q.path_[i]);
	}
	return pnew;
}

