#include "Matrix_Implementation.h"
#include "Vector_Implementation.h"
#include "SingleParticleOperator.h"

typedef complex<double> cd;
typedef double doub;

// Vector instantiations
template class Vector<double>;
template class Vector<complex<double>>;

template Vector<complex<double>> multAB<cd>(const Matrix<complex<double>>& A, const Vector<complex<double>>& B);
template Vector<double> multAB<doub>(const Matrix<double>& A, const Vector<double>& B);
template Vector<cd> multATB<cd>(const Matrix<cd>& A, const Vector<cd>& B);

// Matrix instantiations
template class Matrix<double>;
template class Matrix<complex<double>>;
template class Matrix<bool>;

template Matrix<cd> multAB(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<cd> addAB(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<cd> multscalar<cd, cd>(const cd sca, const Matrix<cd>& B);
template Matrix<cd> multscalar<cd, doub>(const double sca, const Matrix<cd>& B);
template Matrix<cd> substAB<cd>(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<cd> UnitarySimilarityTrafo<cd>(const Matrix<cd>& A, const Matrix<cd>& B);
template Matrix<double> UnitarySimilarityTrafo<doub>(const Matrix<double>& A, const Matrix<double>& B);
template Matrix<cd> Merge(const Matrix<cd>& A, const Matrix<cd>& B,
		const Matrix<cd>& AB);
