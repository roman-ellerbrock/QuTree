//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/FactorMatrix.h"

typedef complex<double> cd;
typedef double d;

// FactorMatrix instantiations
// complex double
template class FactorMatrix<cd>;
template ostream& operator<< <cd>(ostream& os, const FactorMatrix<cd>& A);
template istream& operator>> <cd>(istream& is, FactorMatrix<cd>& A);

// double
template class FactorMatrix<d>;
template ostream& operator<< <d>(ostream& os, const FactorMatrix<d>& A);
template istream& operator>> <d>(istream& is, FactorMatrix<d>& A);
