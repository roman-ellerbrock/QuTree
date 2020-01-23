//
// Created by Roman Ellerbrock on 2020-01-23.
//
#include "TensorTree.h"
#include "TensorTree_Implementation.h"

//template class TensorTree<double>;
typedef complex<double> cd;
template class TensorTree<cd>;

template ostream& operator<< <cd>(ostream& , const TensorTree<cd>& );
template istream& operator>> <cd>(istream& , TensorTree<cd>& );
//template <typename T>
//istream& operator>>(istream& is, TensorTree<T>& t);
