#include "GraphFactory.h"
#include "TensorNetwork.h"
#include "gtest/gtest.h"
#include <ranges>

using namespace qutree;
using namespace std;

TEST(TensorNetwork, shape) {
  Graph graph = balancedBinaryTree(3);
  NetworkShape shape = standardShape(graph, 5, 10);

  ASSERT_EQ(5, shape.nodes_.size());
  ASSERT_EQ(11, shape.edges_.size());

  ASSERT_EQ(std::vector<index_t>({10, 5}), shape.nodes_[0]);
  ASSERT_EQ(std::vector<index_t>({10, 5}), shape.nodes_[1]);
  ASSERT_EQ(std::vector<index_t>({10, 5}), shape.nodes_[2]);
  ASSERT_EQ(std::vector<index_t>({5, 5, 5}), shape.nodes_[3]);
  ASSERT_EQ(std::vector<index_t>({5, 5}), shape.nodes_[4]);

  ASSERT_EQ(std::vector<index_t>({10, 10}), (shape.edges_[{-1, 0}]));
  ASSERT_EQ(std::vector<index_t>({10, 10}), (shape.edges_[{-2, 1}]));
  ASSERT_EQ(std::vector<index_t>({10, 10}), (shape.edges_[{-3, 2}]));

  ASSERT_EQ(std::vector<index_t>({5, 5}), (shape.edges_[{0, 3}]));
  ASSERT_EQ(std::vector<index_t>({5, 5}), (shape.edges_[{1, 3}]));
  ASSERT_EQ(std::vector<index_t>({5, 5}), (shape.edges_[{2, 4}]));
  ASSERT_EQ(std::vector<index_t>({5, 5}), (shape.edges_[{3, 4}]));
}

TEST(TensorNetwork, sizesIntArrayRef) {
  Tensor A = tensorlib::rand({2, 3, 4});
  tensorlib::IntArrayRef shape = A.sizes();
}

// Function that takes a function pointer as an argument
void processTensor(tensorlib::Tensor (*function)(at::IntArrayRef, at::TensorOptions)) {
  tensorlib::IntArrayRef size = {3, 4};
  tensorlib::TensorOptions options = {};
  tensorlib::Tensor tensor = function(size, options);

  // Do something with the generated tensor
  // ...
}

TEST(TensorNetwork, passRand) {
  // Pass torch::rand as the function argument
  processTensor(tensorlib::rand);
  processTensor(tensorlib::zeros);
}

TEST (TensorNetwork, randTN)
{
    Graph graph = balancedBinaryTree(3);
    NetworkShape shape = standardShape(graph, 5, 10);
    TensorNetwork tn = createTN(shape, TN, tensorlib::rand);

    ASSERT_EQ(5, shape.nodes_.size());
    ASSERT_EQ(11, shape.edges_.size());

    ASSERT_EQ(std::vector<index_t>({10, 5}), tn.nodes_[0].sizes());
    ASSERT_EQ(std::vector<index_t>({10, 5}), tn.nodes_[1].sizes());
    ASSERT_EQ(std::vector<index_t>({10, 5}), tn.nodes_[2].sizes());
    ASSERT_EQ(std::vector<index_t>({5, 5, 5}), tn.nodes_[3].sizes());
    ASSERT_EQ(std::vector<index_t>({5, 5}), tn.nodes_[4].sizes());
}

TEST (TensorNetwork, randMT)
{
    Graph graph = balancedBinaryTree(3);
    NetworkShape shape = standardShape(graph, 5, 10);
    TensorNetwork mn = createTN(shape, MN, tensorlib::rand);

    ASSERT_EQ(5, shape.nodes_.size());
    ASSERT_EQ(11, shape.edges_.size());

    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-1, 0}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-2, 1}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-3, 2}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{0, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{1, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{2, 4}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{3, 4}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({0, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({1, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({2, 4})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({3, 4})].sizes()));
}

TEST (TensorNetwork, randCTN)
{
    Graph graph = balancedBinaryTree(3);
    NetworkShape shape = standardShape(graph, 5, 10);
    TensorNetwork ctn = createTN(shape, CTN, tensorlib::rand);

    ASSERT_EQ(5, shape.nodes_.size());
    ASSERT_EQ(11, shape.edges_.size());

    ASSERT_EQ(std::vector<index_t>({10, 5}), ctn.nodes_[0].sizes());
    ASSERT_EQ(std::vector<index_t>({10, 5}), ctn.nodes_[1].sizes());
    ASSERT_EQ(std::vector<index_t>({10, 5}), ctn.nodes_[2].sizes());
    ASSERT_EQ(std::vector<index_t>({5, 5, 5}), ctn.nodes_[3].sizes());
    ASSERT_EQ(std::vector<index_t>({5, 5}), ctn.nodes_[4].sizes());

    ASSERT_EQ(std::vector<index_t>({10, 10}), (ctn.edges_[{-1, 0}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (ctn.edges_[{-2, 1}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (ctn.edges_[{-3, 2}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[{0, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[{1, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[{2, 4}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[{3, 4}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[flip({0, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[flip({1, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[flip({2, 4})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (ctn.edges_[flip({3, 4})].sizes()));
}

TEST (TensorNetwork, randTNfromTN)
{
    Graph graph = balancedBinaryTree(3);
    NetworkShape shape = standardShape(graph, 5, 10);
    TensorNetwork tmp = createTN(shape, TN, tensorlib::rand);
    TensorNetwork tn = createTN(tmp, TN, tensorlib::rand);

    ASSERT_EQ(5, shape.nodes_.size());
    ASSERT_EQ(11, shape.edges_.size());

    ASSERT_EQ(std::vector<index_t>({10, 5}), tn.nodes_[0].sizes());
    ASSERT_EQ(std::vector<index_t>({10, 5}), tn.nodes_[1].sizes());
    ASSERT_EQ(std::vector<index_t>({10, 5}), tn.nodes_[2].sizes());
    ASSERT_EQ(std::vector<index_t>({5, 5, 5}), tn.nodes_[3].sizes());
    ASSERT_EQ(std::vector<index_t>({5, 5}), tn.nodes_[4].sizes());
}

TEST (TensorNetwork, randMTfromMN)
{
    Graph graph = balancedBinaryTree(3);
    NetworkShape shape = standardShape(graph, 5, 10);
    TensorNetwork tmp = createTN(shape, MN, tensorlib::rand);
    TensorNetwork mn = createTN(tmp, MN, tensorlib::rand);

    ASSERT_EQ(5, shape.nodes_.size());
    ASSERT_EQ(11, shape.edges_.size());

    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-1, 0}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-2, 1}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-3, 2}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{0, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{1, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{2, 4}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{3, 4}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({0, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({1, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({2, 4})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({3, 4})].sizes()));
}



TEST (TensorNetwork, randMTfromTN)
{
    Graph graph = balancedBinaryTree(3);
    NetworkShape shape = standardShape(graph, 5, 10);
    TensorNetwork ctn = createTN(shape, CTN, tensorlib::rand);
    TensorNetwork mn = createTN(ctn, MN, tensorlib::rand);

    ASSERT_EQ(5, shape.nodes_.size());
    ASSERT_EQ(11, shape.edges_.size());

    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-1, 0}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-2, 1}].sizes()));
    ASSERT_EQ(std::vector<index_t>({10, 10}), (mn.edges_[{-3, 2}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{0, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{1, 3}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{2, 4}].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[{3, 4}].sizes()));

    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({0, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({1, 3})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({2, 4})].sizes()));
    ASSERT_EQ(std::vector<index_t>({5, 5}), (mn.edges_[flip({3, 4})].sizes()));
}

