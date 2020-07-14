//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Matrix_Implementation.h"
#include "Core/Matrix_Extension_Implementation.h"

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

template Matrix<cd> UnitarySimilarityTrafo<cd>(const Matrix<cd>& A,
	const Matrix<cd>& B);
template Matrix<double> UnitarySimilarityTrafo<doub>(const Matrix<double>& A,
	const Matrix<double>& B);

template Matrix<cd> Re(const Matrix<cd>& A);
template Matrix<d> Re(const Matrix<d>& A);

///////////////////////////////////////////////////////////////
/// Convenience & Management
///////////////////////////////////////////////////////////////

template Matrix<cd> IdentityMatrix(size_t dim);
template Matrix<doub> IdentityMatrix(size_t dim);

template double Residual(const Matrixcd& A, const Matrixcd& B);
template double Residual(const Matrixd& A, const Matrixd& B);

template Matrix<cd> Merge(const Matrix<cd>& A, const Matrix<cd>& B, const Matrix<cd>& AB);
template Matrix<d> Merge(const Matrix<d>& A, const Matrix<d>& B, const Matrix<d>& AB);

template Matrix<cd> Submatrix(const Matrix<cd> A, size_t dim1, size_t dim2);
template Matrix<d> Submatrix(const Matrix<d> A, size_t dim1, size_t dim2);

template ostream& operator<< <cd> (ostream& os, const Matrixcd& A);
template istream& operator>> <cd> (istream& is, Matrixcd& A);

template ostream& operator<< <doub> (ostream& os, const Matrix<doub>& A);
template istream& operator>> <doub> (istream& is, Matrix<doub>& A);

template Matrix<double> EuclideanDistance(const Matrix<double>& A);

///////////////////////////////////////////////////////////////
/// Matrix Extension
///////////////////////////////////////////////////////////////

template SpectralDecompositioncd sqrt(SpectralDecompositioncd X);
template SpectralDecompositiond sqrt(SpectralDecompositiond X);

template SpectralDecompositioncd inverse(SpectralDecompositioncd X, double eps);
template SpectralDecompositiond inverse(SpectralDecompositiond X, double eps);

template Matrixcd toMatrix(const SpectralDecompositioncd& X);
template Matrixd toMatrix(const SpectralDecompositiond& X);

template Matrixcd BuildInverse(const SpectralDecompositioncd& X, double eps);
template Matrixd BuildInverse(const SpectralDecompositiond& X, double eps);

template SpectralDecomposition<double> reduceRank(const SpectralDecomposition<double>& x, size_t rank);
template SpectralDecomposition<cd> reduceRank(const SpectralDecomposition<cd>& x, size_t rank);

