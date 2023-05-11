#include "gtest/gtest.h"
#include "GraphFactory.h"
#include "TensorNetwork.h"
#include <ranges>

using namespace qutree;
using namespace std;

TEST (TensorNetwork, shape)
{
    Graph graph = balancedBinaryTree(3);
    NetworkShape shape = standardShape(graph, 5, 10);

    cout << graph << endl;

    ASSERT_EQ(5, shape.nodes_.size());
    ASSERT_EQ(8, shape.edges_.size());
    ASSERT_EQ(3, shape.leaves_.size());

    // bottom layer
    // have to add leave edges to bottomlayer
//    ASSERT_EQ(std::vector<index_t>({10, 5}), shape.nodes_[0]);
//    ASSERT_EQ(std::vector<index_t>({10, 5}), shape.nodes_[1]);
//    ASSERT_EQ(std::vector<index_t>({10, 5}), shape.nodes_[2]);
    // upper
    ASSERT_EQ(std::vector<index_t>({5, 5, 5}), shape.nodes_[3]);
    ASSERT_EQ(std::vector<index_t>({5, 5}), shape.nodes_[4]);

    cout << "Nodes:\n";
    for (pair<Node, vector<index_t>> p : shape.nodes_) {
        cout << p.first << "\t";
        for (auto x : p.second) {
            cout << x << " ";
        } cout << endl;
    }

    cout << "Edges:\n";
    for (pair<Edge, vector<index_t>> p : shape.edges_) {
        cout << p.first << "\t";
        for (auto x : p.second) {
            cout << x << " ";
        } cout << endl;
    }

    cout << "Leaves:\n";
    for (pair<Edge, vector<index_t>> p : shape.leaves_) {
        cout << p.first << "\t";
        for (auto x : p.second) {
            cout << x << " ";
        } cout << endl;
    }

}

TEST (TensorNetwork, construct)
{
}
