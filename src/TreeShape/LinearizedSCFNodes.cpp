#include "TreeShape/LinearizedSCFNodes.h"

void LinearizedSCFNodes::linearize(Tree& tree) {
    AbstractNode* next_node = tree.topNode().nextSCFNode(nullptr);
    while(next_node != nullptr){
        pointers_.push_back(reinterpret_cast<Node*>(next_node));
        addresses_.push_back(reinterpret_cast<Node*>(next_node)->address());
        next_node = tree.topNode().nextSCFNode(next_node);
    }

    tree.topNode().resetSCFstatus();
}

const std::vector<int>& LinearizedSCFNodes::getAddresses() const {
    return addresses_;
}

const std::vector<Node *>& LinearizedSCFNodes::getNodes() const {
    return pointers_;
}

void LinearizedSCFNodes::clear() {
    addresses_.clear();
    pointers_.clear();
}
