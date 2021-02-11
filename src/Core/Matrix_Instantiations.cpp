//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Matrix_Implementation.h"
#include "Core/Matrix_Extension_Implementation.h"
#include "Util/RandomProjector_Implementation.h"

typedef complex<double> cd;
typedef double doub;
typedef double d;

///////////////////////////////////////////////////////////////
/// Matrix class instantiations
///////////////////////////////////////////////////////////////

template class Matrix<int>;
template class Matrix<double>;
template class Matrix<complex<double>>;

///////////////////////////////////////////////////////////////
/// Arithmetic
///////////////////////////////////////////////////////////////

template Vector<complex<double>> multAB<cd>(const Matrix<complex<double>>& A, const Vector<complex<double>>& B);
template Vector<double> multAB<doub>(const Matrix<double>& A, const Vector<double>& B);

template Vector<cd> multATB<cd>(const Matrix<cd>& A, const Vector<cd>& B);
template Vector<d> multATB<d>(const Matrix<d>& A, const Vector<d>& B);

template Matrix<double> multAB(const Matrix<double>& A, const Matrix<double>& B);
template Matrix<cd> multAB(const Matrix<cd>& A, const Matrix<cd>& B);

template Matrix<double> multATB(const Matrix<double>& A, const Matrix<double>& B);
template Matrix<cd> multATB(const Matrix<cd>& A, const Matrix<cd>& B);

template Matrix<cd> addAB(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<doub> addAB(const Matrix<doub>& A, const Matrix<doub>& B);

template Matrix<doub> substAB(const Matrix<doub>& A, const Matrix<doub>& B);
template Matrix<cd> substAB<cd>(const Matrix<cd>& A, const Matrix<cd>& B);

template Matrix<cd> multscalar<cd, cd>(const cd sca, const Matrix<cd>& B);
template Matrix<cd> multscalar<cd, doub>(const double sca, const Matrix<cd>& B);
template Matrix<doub> multscalar<doub, doub>(const double sca,
	const Matrix<doub>& B);

template Matrix<cd> unitarySimilarityTrafo<cd>(const Matrix<cd>& A,
	const Matrix<cd>& B);
template Matrix<double> unitarySimilarityTrafo<doub>(const Matrix<double>& A,
	const Matrix<double>& B);

template Matrix<cd> re(const Matrix<cd>& A);
template Matrix<d> re(const Matrix<d>& A);

///////////////////////////////////////////////////////////////
/// Convenience & Management
///////////////////////////////////////////////////////////////

template Matrix<cd> identityMatrix(size_t dim);
template Matrix<doub> identityMatrix(size_t dim);

template double residual(const Matrixcd& A, const Matrixcd& B);
template double residual(const Matrixd& A, const Matrixd& B);

template Matrix<cd> Regularize(const Matrix<cd>& A, double eps);
template Matrix<d> Regularize(const Matrix<d>& A, double eps);

template Matrix<cd> merge(const Matrix<cd>& A, const Matrix<cd>& B, const Matrix<cd>& AB);
template Matrix<d> merge(const Matrix<d>& A, const Matrix<d>& B, const Matrix<d>& AB);

template Matrix<cd> subMatrix(const Matrix<cd> A, size_t dim1, size_t dim2);
template Matrix<d> subMatrix(const Matrix<d> A, size_t dim1, size_t dim2);

template ostream& operator<< <cd> (ostream& os, const Matrixcd& A);
template istream& operator>> <cd> (istream& is, Matrixcd& A);

template ostream& operator<< <doub> (ostream& os, const Matrix<doub>& A);
template istream& operator>> <doub> (istream& is, Matrix<doub>& A);

template Matrix<double> euclideanDistance(const Matrix<double>& A);

///////////////////////////////////////////////////////////////
/// Matrix Extension
///////////////////////////////////////////////////////////////

template SpectralDecompositioncd sqrt(SpectralDecompositioncd X);
template SpectralDecompositiond sqrt(SpectralDecompositiond X);

template SpectralDecompositioncd inverse(SpectralDecompositioncd X, double eps);
template SpectralDecompositiond inverse(SpectralDecompositiond X, double eps);

template Matrixcd toMatrix(const SpectralDecompositioncd& X);
template Matrixd toMatrix(const SpectralDecompositiond& X);

//template Matrixcd inverse(const SpectralDecompositioncd& X, double eps);
//template Matrixd inverse(const SpectralDecompositiond& X, double eps);

template SpectralDecomposition<double> reduceRank(const SpectralDecomposition<double>& x, size_t rank);
template SpectralDecomposition<cd> reduceRank(const SpectralDecomposition<cd>& x, size_t rank);

template SpectralDecomposition<cd>
    Random::DiagonalizeRandom<cd, Matrixcd>(const Matrixcd& A,
	size_t rank, size_t pow, mt19937& gen);
