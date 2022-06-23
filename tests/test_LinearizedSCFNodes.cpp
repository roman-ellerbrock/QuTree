//
// Created by Hannes Hoppe on 22.06.22.
//

#include <gtest/gtest.h>
#include "TreeShape/LinearizedSCFNodes.h"
#include "TreeShape/TreeFactory.h"

TEST(LinearizedSCFNode, correct_order){

    Tree tree = TreeFactory::balancedTree(543, 4, 11);
    LinearizedSCFNodes test;
    test.linearize(tree);

    // prepare list to compare against
    std::vector<Node*> node_list;
    std::vector<int> address_list;

    AbstractNode* nextnode = tree.topNode().nextSCFNode(nullptr);
    while(nextnode != nullptr){
        node_list.push_back(reinterpret_cast<Node*>(nextnode));
        address_list.push_back(reinterpret_cast<Node*>(nextnode)->address());
        nextnode = tree.topNode().nextSCFNode(nextnode);
    }

    const std::vector<Node*> test_node_list = test.getNodes();
    const std::vector<int> test_address_list = test.getAddresses();

    ASSERT_EQ(test_address_list.size(), address_list.size());
    ASSERT_EQ(test_node_list.size(), node_list.size());
    ASSERT_EQ(test_node_list.size(), test_address_list.size());

    for(size_t i = 0; i < node_list.size(); ++i){
        ASSERT_EQ(test_node_list[i], node_list[i]);
        ASSERT_EQ(test_address_list[i], address_list[i]);
    }
}

TEST(LinearizedSCFNode, clear){
    Tree tree = TreeFactory::balancedTree(543, 4, 11);
    LinearizedSCFNodes test;
    test.linearize(tree);

    ASSERT_EQ(test.getAddresses().size(), 2 * tree.nodes().size() - 1);
    ASSERT_EQ(test.getNodes().size(), 2 * tree.nodes().size() - 1);

    test.clear();
    ASSERT_EQ(test.getAddresses().size(), 0);
    ASSERT_EQ(test.getNodes().size(), 0);

}