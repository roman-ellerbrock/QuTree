//
// Created by Hannes Hoppe on 17.06.22.
//

#include <gtest/gtest.h>
#include "TreeClasses/MatrixTensorTree.h"
#include "TreeShape/TreeFactory.h"
#include "TreeShape/Node.h"

TEST(Node, correctSCFNodeOrder){

    Tree tree = TreeFactory::balancedTree(7, 2, 4);

    // predicted order
    const std::vector<int> predicted_order{9,3,0,3,1,3,2,3,9,8,4,8,5,8,6,8,7,8,9};

    std::vector<AbstractNode*> node_list;
    AbstractNode* nextnode = tree.topNode().nextSCFNode(nullptr);

    while(nextnode != nullptr){
        node_list.push_back(nextnode);
        nextnode = tree.topNode().nextSCFNode(nextnode);
    }

    ASSERT_EQ(predicted_order.size(), node_list.size());
    for(size_t i = 0; i < node_list.size(); ++i){
        ASSERT_EQ(predicted_order[i], reinterpret_cast<Node*>(node_list[i])->address());
    }
}

TEST(Node, nextSCFnodeWorksMultipleTimes){

    Tree tree = TreeFactory::balancedTree(543, 4, 11);

    std::vector<int> nodes_visited;
    AbstractNode* next_node = tree.topNode().nextSCFNode(nullptr);

    // build up node list to compare against
    while(next_node != nullptr){
        nodes_visited.push_back(reinterpret_cast<Node*>(next_node)->address());
        next_node = tree.topNode().nextSCFNode(next_node);
    }

    ASSERT_EQ(nodes_visited.size(), 2 * tree.nodes().size() - 1);

    AbstractNode* prev_node = nullptr;
    for(const auto& i : nodes_visited){
        prev_node = tree.topNode().nextSCFNode(prev_node);
        ASSERT_EQ(i, reinterpret_cast<Node*>(prev_node)->address());
    }
    ASSERT_EQ(tree.topNode().nextSCFNode(prev_node), nullptr);

    prev_node = nullptr;
    for(const auto& i : nodes_visited){
        prev_node = tree.topNode().nextSCFNode(prev_node);
        ASSERT_EQ(i, reinterpret_cast<Node*>(prev_node)->address());
    }
    ASSERT_EQ(tree.topNode().nextSCFNode(prev_node), nullptr);
}