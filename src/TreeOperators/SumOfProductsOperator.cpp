#include "TreeOperators/SumOfProductsOperator_Implementation.h"
#include "TreeOperators/SOPVector.h"


typedef complex<double> cd;
template class SumOfProductsOperator<cd>;

template SumOfProductsOperator<cd> operator*<cd>(cd c, const SOP<cd>& A);
template SumOfProductsOperator<cd> operator*(const SOP<cd>& A, complex<double> c);

template SumOfProductsOperator<cd> operator*(const MLO<cd>& M, const SOP<cd>& A);
template SumOfProductsOperator<cd> operator*(const SOP<cd>& A, const MLO<cd>& M);
template SumOfProductsOperator<cd> operator*(const SOP<cd>& A, const SOP<cd>& B);

template SumOfProductsOperator<cd> operator+(const SOP<cd>& A, const SOP<cd>& B);

/// double
typedef double d;
template class SumOfProductsOperator<d>;
template SumOfProductsOperator<d> operator*<d>(d c, const SOP<d>& A);
template SumOfProductsOperator<d> operator*(const SOP<d>& A, double c);

template SumOfProductsOperator<d> operator*(const MLO<d>& M, const SOP<d>& A);
template SumOfProductsOperator<d> operator*(const SOP<d>& A, const MLO<d>& M);
template SumOfProductsOperator<d> operator*(const SOP<d>& A, const SOP<d>& B);

template SumOfProductsOperator<d> operator+(const SOP<d>& A, const SOP<d>& B);
