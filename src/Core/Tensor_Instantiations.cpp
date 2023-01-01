//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Tensor_Implementation.h"

typedef complex<double> cd;
typedef double doub;
typedef double d;

// Tensor instantiations
template class Tensor<double>;
template class Tensor<complex<double>>;

template Tensor<cd> productElementwise(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<d> productElementwise(const Tensor<d>& A, const Tensor<d>& B);

//template Matrix<cd> HoleProduct<cd>(const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
//template void HoleProduct<cd>(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
template Tensor<cd> matrixTensor<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template void matrixTensor<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Tensor<cd> tMatrixTensor<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> multStateAB(const Matrix<cd>& A, const Tensor<cd>& B);
template Tensor<cd> multStateArTB<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B);
template void multStateArTB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B);
template void gramSchmidt<cd>(Tensor<cd>& A);
template void multStateAB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, bool zero);
template Tensor<cd> project<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> projectOut<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> projectOrthogonal<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template void multAdd<cd, cd>(Tensor<cd>& A, const Tensor<cd>& B, cd coeff);
template Tensor<cd> conj<cd>(Tensor<cd> A);
template double residual(Tensorcd A, const Tensorcd& B);
template Matrix<cd> toMatrix(const Tensor<cd>& A);
template Matrix<cd> moveToMatrix(Tensor<cd>& A);
template Tensor<cd> toTensor(const Matrix<cd>& B);
template Matrix<cd> toMatrix(const Tensor<cd>& A, size_t mode);
template Tensor<cd> toTensor(const Matrix<cd>& B, const TensorShape& shape, size_t mode);
template Tensor<cd> moveToTensor(Matrix<cd>& B);
template void contraction<cd>(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t before, size_t active1, size_t active2, size_t behind);
template void matrixTensor<cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>&  B, size_t before, size_t activeC, size_t activeB, size_t after, bool zero);
template Matrix<cd> contraction(const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
template void contraction(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t k, bool zero);
template Tensor<cd> qr(const Tensor<cd>& A);
template Tensor<cd> qr(const Tensor<cd>& A, size_t mode);

template void tensorMatrix(Tensor<cd>& C, const Tensor<cd>& B, const Matrix<cd>& A, size_t mode, bool zero);
template Tensor<cd> tensorMatrix(const Tensor<cd>& B, const Matrix<cd>& A, size_t mode);

template ostream& operator<< <cd> (ostream&, const Tensor<cd>& );
template istream& operator>> <cd> (istream&, Tensor<cd>& );
template bool operator== <cd>(const Tensor<cd>& A, const Tensor<cd>& B);

//template Matrix<double> HoleProduct<double>(const Tensor<double>& A, const Tensor<double>& B, size_t k);
//template void HoleProduct<doub>(Matrix<doub>& S, const Tensor<doub>& A, const Tensor<doub>& B, size_t k);
template Tensor<double> matrixTensor<doub, doub>(const Matrix<double>& A, const Tensor<double>& B, size_t mode);
template void matrixTensor<doub, doub>(Tensor<double>& C, const Matrix<double>& A, const Tensor<double>& B, size_t mode, bool zero);
template void tensorMatrix(Tensor<double>& C, const Tensor<double>& B, const Matrix<double>& A, size_t mode, bool zero);
template Tensor<double> tensorMatrix(const Tensor<double>& B, const Matrix<double>& A, size_t mode);
template Tensor<double> tMatrixTensor<doub, doub>(const Matrix<double>& A, const Tensor<double>& B, size_t mode);
template Tensor<double> multStateAB(const Matrix<double>& A, const Tensor<double>& B);
template Tensor<d> projectOut<d>(const Tensor<d>& A, const Tensor<d>& B);
template Tensor<double> multStateArTB<doub, doub>(const Matrix<double>& A, const Tensor<double>& B);
template void multStateArTB<doub, doub>(Tensor<doub>& C, const Matrix<doub>& A, const Tensor<doub>& B);
template void gramSchmidt<double>(Tensor<double>& A);
template void multStateAB<doub, doub>(Tensor<double>& C, const Matrix<double>& A, const Tensor<double>& B, bool zero);
template void multAdd<doub, doub>(Tensor<double>& A, const Tensor<double>& B, doub coeff);
template Tensor<double> conj<double>(Tensor<double> A);
template double residual(Tensord A, const Tensord& B);
template Matrix<double> toMatrix(const Tensor<double>& A);
template Matrix<d> moveToMatrix(Tensor<d>& A);
template Tensor<double> toTensor(const Matrix<double>& B);
template void contraction<doub>(Matrix<doub>& S, const Tensor<doub>& A, const Tensor<doub>& B, size_t before, size_t active1, size_t active2, size_t behind);
template void matrixTensor<doub>(Tensor<doub>& C, const Matrix<doub>& A, const Tensor<doub>&  B, size_t before, size_t activeC, size_t activeB, size_t after, bool zero);
template Matrix<doub> contraction(const Tensor<doub>& A, const Tensor<doub>& B, size_t k);
template void contraction(Matrix<doub>& S, const Tensor<doub>& A, const Tensor<doub>& B, size_t k, bool zero);
template Tensor<doub> qr(const Tensor<doub>& A);
template Tensor<doub> qr(const Tensor<doub>& A, size_t mode);

template ostream& operator<< <doub> (ostream&, const Tensor<doub>& );
template istream& operator>> <doub> (istream&, Tensor<doub>& );
template bool operator== <doub>(const Tensor<doub>& A, const Tensor<doub>& B);

/// Mixed double/complex<double>
template void multAdd<cd, double>(Tensor<cd>& A, const Tensor<cd>& B, double coeff);
template Tensor<cd> matrixTensor<cd, doub>(const Matrix<doub>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> tMatrixTensor<cd, doub>(const Matrix<doub>& A, const Tensor<cd>& B, size_t mode);

template Tensord ones<double>(const TensorShape& shape);
template Tensorcd ones<complex<double>>(const TensorShape& shape);

template Tensord rand<double>(const TensorShape& shape);
template Tensorcd rand<complex<double>>(const TensorShape& shape);
