//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/FactorMatrix.h"

typedef complex<double> cd;
typedef double d;

// FactorMatrix instantiations
// complex double
template class FactorMatrix<cd>;
template void HoleProduct(FactorMatrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t mode);
template FactorMatrix<cd> HoleProduct(const Tensor<cd>& A, const Tensor<cd>& B, size_t mode);
template void DotProduct(FactorMatrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B);
template ostream& operator<< <cd>(ostream& os, const FactorMatrix<cd>& A);
template istream& operator>> <cd>(istream& is, FactorMatrix<cd>& A);

// double
template class FactorMatrix<d>;
template void HoleProduct(FactorMatrix<d>& S, const Tensor<d>& A, const Tensor<d>& B, size_t mode);
template FactorMatrix<d> HoleProduct(const Tensor<d>& A, const Tensor<d>& B, size_t mode);
template void DotProduct(FactorMatrix<d>& S, const Tensor<d>& A, const Tensor<d>& B);
template ostream& operator<< <d>(ostream& os, const FactorMatrix<d>& A);
template istream& operator>> <d>(istream& is, FactorMatrix<d>& A);

