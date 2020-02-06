//
// Created by Roman Ellerbrock on 2020-01-23.
//
#include "Tree/TensorTree.h"
#include "Tree/TensorTree_Implementation.h"

//template class Tree<double>;
typedef complex<double> cd;
template class TensorTree<cd>;

template ostream& operator<< <cd>(ostream& , const TensorTree<cd>& );
template istream& operator>> <cd>(istream& , TensorTree<cd>& );

typedef double d;
template class TensorTree<d>;

template ostream& operator<< <d>(ostream& , const TensorTree<d>& );
template istream& operator>> <d>(istream& , TensorTree<d>& );
