//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Tensor_Implementation.h"

typedef complex<double> cd;
typedef double doub;

// Tensor instantiations
template class Tensor<double>;
template class Tensor<complex<double>>;

template Matrix<cd> HoleProduct<cd>(const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
template Tensor<cd> multAB<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> multAB<cd, doub>(const Matrix<doub>& A, const Tensor<cd>& B, size_t mode);
template void multAB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Tensor<cd> multATB<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> multATB<cd, doub>(const Matrix<doub>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> multStateAB(const Matrix<cd>& A, const Tensor<cd>& B);
template Tensor<cd> multStateArTB<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B);
template void GramSchmidt<cd>(Tensor<cd>& A);
template void multStateAB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, bool zero);
template Tensor<cd> Project<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> ProjectOut<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> ProjectOrthogonal<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template void multAdd<cd, cd>(Tensor<cd>& A, const Tensor<cd>& B, cd coeff);
template void multAdd<cd, double>(Tensor<cd>& A, const Tensor<cd>& B, double coeff);
template Tensor<cd> conj<cd>(Tensor<cd> A);
template double Residual(Tensorcd A, const Tensorcd& B);
template double Residual(Tensord A, const Tensord& B);
