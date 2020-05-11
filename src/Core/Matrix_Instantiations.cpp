//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Matrix_Implementation.h"

typedef complex<double> cd;

// Matrix instantiations
template class Matrix<double>;
template class Matrix<complex<double>>;

template Vector<complex<double>> multAB<cd>(const Matrix<complex<double>>& A, const Vector<complex<double>>& B);
template Vector<cd> multATB<cd>(const Matrix<cd>& A, const Vector<cd>& B);

template Matrix<cd> multAB(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<cd> addAB(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<cd> multscalar<cd, cd>(const cd sca, const Matrix<cd>& B);
template Matrix<cd> substAB<cd>(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<cd> IdentityMatrix(size_t dim);
template Matrix<cd> UnitarySimilarityTrafo<cd>(const Matrix<cd>& A,
	const Matrix<cd>& B);
template Matrix<cd> Merge(const Matrix<cd>& A, const Matrix<cd>& B,
	const Matrix<cd>& AB);
template double Residual(const Matrixcd& A, const Matrixcd& B);
template double Residual(const Matrixd& A, const Matrixd& B);
template Matrix<cd> Submatrix(const Matrix<cd> A, size_t dim1, size_t dim2);


template ostream& operator<< <cd> (ostream& os, const Matrixcd& A);
template istream& operator>> <cd> (istream& is, Matrixcd& A);

typedef double doub;
template Vector<double> multAB<doub>(const Matrix<double>& A,
	const Vector<double>& B);
template Matrix<double> multAB(const Matrix<double>& A, const Matrix<double>& B);
template Matrix<doub> addAB(const Matrix<doub>& A, const Matrix<doub>& B);
template Matrix<doub> substAB(const Matrix<doub>& A, const Matrix<doub>& B);
template Matrix<cd> multscalar<cd, doub>(const double sca,
	const Matrix<cd>& B);
template Matrix<doub> multscalar<doub, doub>(const double sca,
	const Matrix<doub>& B);
template Matrix<doub> IdentityMatrix(size_t dim);
template Matrix<double> UnitarySimilarityTrafo<doub>(const Matrix<double>& A,
	const Matrix<double>& B);
template ostream& operator<< <doub> (ostream& os, const Matrix<doub>& A);
template istream& operator>> <doub> (istream& is, Matrix<doub>& A);


template Matrix<double> EuclideanDistance(const Matrix<double>& A);

template SpectralDecompositioncd sqrt(SpectralDecompositioncd X);
template SpectralDecompositiond sqrt(SpectralDecompositiond X);
template SpectralDecompositioncd inverse(SpectralDecompositioncd X, double eps);
template SpectralDecompositiond inverse(SpectralDecompositiond X, double eps);
template Matrixcd BuildMatrix(const SpectralDecompositioncd& X);
template Matrixd BuildMatrix(const SpectralDecompositiond& X);
template Matrixcd BuildInverse(const SpectralDecompositioncd& X, double eps);
template Matrixd BuildInverse(const SpectralDecompositiond& X, double eps);
template Matrix<doub> Submatrix(const Matrix<doub> A, size_t dim1, size_t dim2);

