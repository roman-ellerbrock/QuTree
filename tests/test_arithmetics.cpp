#include "gtest/gtest.h"
#include "TensorNetwork.h"
#include "GraphFactory.h"
#include "arithmetics.h"


using namespace qutree;
using namespace std;

TEST(Arithmetics, dot)
{
  Graph graph = balancedBinaryTree(4);
  NetworkShape shape = standardShape(graph, 3, 5);
  TensorNetwork tn = createTN(shape, CTN);
  TensorNetwork mt = dotProduct(tn, tn);
  cout << mt << endl;
}
