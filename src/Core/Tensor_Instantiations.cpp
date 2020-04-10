//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Tensor_Implementation.h"

typedef complex<double> cd;

// Tensor instantiations
template class Tensor<double>;
template class Tensor<complex<double>>;

//template Matrix<cd> HoleProduct<cd>(const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
//template void HoleProduct<cd>(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
template Tensor<cd> MatrixTensor<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template void MatrixTensor<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Tensor<cd> multATB<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> multStateAB(const Matrix<cd>& A, const Tensor<cd>& B);
template Tensor<cd> multStateArTB<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B);
template void multStateArTB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B);
template void GramSchmidt<cd>(Tensor<cd>& A);
template void multStateAB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, bool zero);
template Tensor<cd> Project<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> ProjectOut<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> ProjectOrthogonal<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template void multAdd<cd, cd>(Tensor<cd>& A, const Tensor<cd>& B, cd coeff);
template Tensor<cd> conj<cd>(Tensor<cd> A);
template double Residual(Tensorcd A, const Tensorcd& B);
template void TensorContraction<cd>(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t before, size_t active1, size_t active2, size_t behind);
template void MatrixTensor<cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>&  B, size_t before, size_t activeC, size_t activeB, size_t after, bool zero);
template Matrix<cd> Contraction(const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
template void Contraction(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t k, bool zero);

template void TensorMatrix(Tensor<cd>& C, const Tensor<cd>& B, const Matrix<cd>& A, size_t mode, bool zero);
template Tensor<cd> TensorMatrix(const Tensor<cd>& B, const Matrix<cd>& A, size_t mode);

template ostream& operator<< <cd> (ostream&, const Tensor<cd>& );
template istream& operator>> <cd> (istream&, Tensor<cd>& );
template bool operator== <cd>(const Tensor<cd>& A, const Tensor<cd>& B);

typedef double doub;
//template Matrix<double> HoleProduct<double>(const Tensor<double>& A, const Tensor<double>& B, size_t k);
//template void HoleProduct<doub>(Matrix<doub>& S, const Tensor<doub>& A, const Tensor<doub>& B, size_t k);
template Tensor<double> MatrixTensor<doub, doub>(const Matrix<double>& A, const Tensor<double>& B, size_t mode);
template void MatrixTensor<doub, doub>(Tensor<double>& C, const Matrix<double>& A, const Tensor<double>& B, size_t mode, bool zero);
template void TensorMatrix(Tensor<double>& C, const Tensor<double>& B, const Matrix<double>& A, size_t mode, bool zero);
template Tensor<double> TensorMatrix(const Tensor<double>& B, const Matrix<double>& A, size_t mode);
template Tensor<double> multATB<doub, doub>(const Matrix<double>& A, const Tensor<double>& B, size_t mode);
template Tensor<double> multStateAB(const Matrix<double>& A, const Tensor<double>& B);
template Tensor<double> multStateArTB<doub, doub>(const Matrix<double>& A, const Tensor<double>& B);
template void multStateArTB<doub, doub>(Tensor<doub>& C, const Matrix<doub>& A, const Tensor<doub>& B);
template void GramSchmidt<double>(Tensor<double>& A);
template void multStateAB<doub, doub>(Tensor<double>& C, const Matrix<double>& A, const Tensor<double>& B, bool zero);
template void multAdd<doub, doub>(Tensor<double>& A, const Tensor<double>& B, doub coeff);
template Tensor<double> conj<double>(Tensor<double> A);
template double Residual(Tensord A, const Tensord& B);
template void TensorContraction<doub>(Matrix<doub>& S, const Tensor<doub>& A, const Tensor<doub>& B, size_t before, size_t active1, size_t active2, size_t behind);
template void MatrixTensor<doub>(Tensor<doub>& C, const Matrix<doub>& A, const Tensor<doub>&  B, size_t before, size_t activeC, size_t activeB, size_t after, bool zero);
template Matrix<doub> Contraction(const Tensor<doub>& A, const Tensor<doub>& B, size_t k);
template void Contraction(Matrix<doub>& S, const Tensor<doub>& A, const Tensor<doub>& B, size_t k, bool zero);

template ostream& operator<< <doub> (ostream&, const Tensor<doub>& );
template istream& operator>> <doub> (istream&, Tensor<doub>& );
template bool operator== <doub>(const Tensor<doub>& A, const Tensor<doub>& B);

/// Mixed double/complex<double>
template void multAdd<cd, double>(Tensor<cd>& A, const Tensor<cd>& B, double coeff);
template Tensor<cd> MatrixTensor<cd, doub>(const Matrix<doub>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> multATB<cd, doub>(const Matrix<doub>& A, const Tensor<cd>& B, size_t mode);


