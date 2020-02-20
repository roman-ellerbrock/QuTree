#include "TreeOperators/SumOfProductsOperator_Implementation.h"


typedef complex<double> cd;
template class SumOfProductsOperator<cd>;

template SumOfProductsOperator<cd> operator*<cd>(cd c, const SOP<cd>& A);
template SumOfProductsOperator<cd> operator*(const SOP<cd>& A, complex<double> c);

template SumOfProductsOperator<cd> operator*(const MLO<cd>& M, const SOP<cd>& A);
template SumOfProductsOperator<cd> operator*(const SOP<cd>& A, const MLO<cd>& M);

template SumOfProductsOperator<cd> operator+(const SOP<cd>& A, const SOP<cd>& B);
