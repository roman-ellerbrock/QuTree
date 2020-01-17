#include "Tensor_Implementation.h"
#include "Matrix_Implementation.h"
#include "Vector_Implementation.h"
#include "SingleParticleOperator.h"
#include "Tensor_Extension_Implementation.h"

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

// Tensor-Extension instantiations
template Matrix<cd> Tensor_Extension::OuterProduct(const Tensor<cd>& A, const Tensor<cd>& B);
template void Tensor_Extension::OuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B);
template Matrix<cd> Tensor_Extension::WeightedOuterProduct(const Tensor<cd>& A, const Tensor<cd>& B,
		const Matrix<cd>& m);
template void Tensor_Extension::WeightedOuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B,
		const Matrix<cd>& rho);
template Tensor<cd> Tensor_Extension::Merge(Tensor<cd> A, const Tensor<cd>& B);
template Matrix<cd> Tensor_Extension::Map(const Tensor<cd>& A);

// SingleParticleOperator instantiations
template class SingleParticleOperator<complex<double>>;
template class SingleParticleOperator<double>;

