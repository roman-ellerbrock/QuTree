//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Vector_Implementation.h"

typedef complex<double> cd;
typedef double d;

// Vector instantiations
template class Vector<int>;
template class Vector<double>;
template class Vector<complex<double>>;

template void normalize(Vectord& a);
template void normalize(Vectorcd& a);

template double residual(const Vectord& A, const Vectord& B);
template double residual(const Vectorcd& A, const Vectorcd& B);

template Vector<cd> reverse(const Vector<cd>& A);
template Vector<d> reverse(const Vector<d>& A);

template cd sum(const Vector<cd>& A);
template d sum(const Vector<d>& A);

template Vector<d> regularize(Vector<d> A, double eps);
template Vector<cd> regularize(Vector<cd> A, double eps);

template Vector<d> inverse(Vector<d> A, d eps);
template Vector<cd> inverse(Vector<cd> A, d eps);
