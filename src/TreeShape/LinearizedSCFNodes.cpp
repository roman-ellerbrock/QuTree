#include "TreeShape/LinearizedSCFNodes.h"

void LinearizedSCFNodes::linearize(Tree& tree) {
    AbstractNode* next_node = tree.topNode().nextSCFNode(nullptr);
    while(next_node != nullptr){
        pointers_.push_back(reinterpret_cast<Node*>(next_node));
        addresses_.push_back(reinterpret_cast<Node*>(next_node)->address());
        next_node = tree.topNode().nextSCFNode(next_node);
    }

    // call to "reset" tree structure into previous state
    tree.topNode().resetSCFstatus();
}

std::vector<int> LinearizedSCFNodes::getAddresses() const {
    return addresses_;
}

std::vector<Node *> LinearizedSCFNodes::getNodes() const {
    return pointers_;
}
