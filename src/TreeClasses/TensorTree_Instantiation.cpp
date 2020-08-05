//
// Created by Roman Ellerbrock on 2020-01-23.
//
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/TensorTree_Implementation.h"
#include "TreeClasses/TensorTreeFunctions_Implementation.h"

//template class Tree<double>;
typedef complex<double> cd;
template class TensorTree<cd>;

template ostream& operator<< <cd>(ostream& , const TensorTree<cd>& );
template istream& operator>> <cd>(istream& , TensorTree<cd>& );

template TensorTree<cd> operator*(cd c, TensorTree<cd> R);
template TensorTree<cd> operator/(cd c, TensorTree<cd> R);

typedef double d;
template class TensorTree<d>;

template ostream& operator<< <d>(ostream& , const TensorTree<d>& );
template istream& operator>> <d>(istream& , TensorTree<d>& );

template TensorTree<d> operator*(d c, TensorTree<d> R);
template TensorTree<d> operator/(d c, TensorTree<d> R);

template void Orthogonal<cd>(TensorTree<cd>& Psi, const Tree& tree);
template void Orthogonal<d>(TensorTree<d>& Psi, const Tree& tree);

template void Orthonormal<cd>(TensorTree<cd>& Psi, const Tree& tree);
template void Orthonormal<d>(TensorTree<d>& Psi, const Tree& tree);

template void TreeFunctions::Adjust(TensorTree<cd>& Psi, Tree& tree, const SpectralDecompositionTree<cd>& X, double eps);
template void TreeFunctions::Adjust(TensorTree<d>& Psi, Tree& tree, const SpectralDecompositionTree<d>& X, double eps);

template void TreeFunctions::Adjust<cd>(TensorTree<cd>& Psi, const Tree& newtree);
template void TreeFunctions::Adjust<d>(TensorTree<d>& Psi, const Tree& newtree);

template void TreeFunctions::Sum(TensorTree<d>& Psi, Tree& tree, const TensorTree<d>& Chi, bool sameLeafs, bool sumToplayer);
template void TreeFunctions::Sum(TensorTree<cd>& Psi, Tree& tree, const TensorTree<cd>& Chi, bool sameLeafs, bool sumToplayer);
