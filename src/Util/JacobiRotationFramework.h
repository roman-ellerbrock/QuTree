#pragma once
#include "Tensor/Tensor.h"

namespace JacobiRotationFramework
/**
 * \namespace JacobiRotationFramework
 * \ingroup Util
 * \brief This namespace contains functions for Jacobi rotations algorithms
 */
{
	/**
	 * \brief Perform a Givens-rotation on an FactorMatrix
	 * @param B matrix to rotate
	 * @param c givins matrix cos element
	 * @param s givins matrix sin element
	 * @param i target index in B
	 * @param j target index in B
	 */
	void givensRotation(Matrixcd& B, complex<double> c, complex<double> s,
		int i, int j);

	/**
	 * \brief Perform a Givens-Rotation on a set of FactorMatrices A
	 * @param A Set of rotated matrices
	 * @param c Givens matrix cos element
	 * @param s Givens matrix sin element
	 * @param i Target index in A
	 * @param j Target index in A
	 */
	void rotateMatrices(vector<Matrixcd>& A, complex<double> c,
		complex<double> s, int i, int j);

	/**
	 * \brief rotate Transformation matrix
	 * @param trafo transformation matrix
	 * @param c Givens matrix cos element
	 * @param s Givens matrix sin element
	 * @param i Target index in A
	 * @param j Target index in A
	 */
	void givensTrafoRotation(Matrixcd& trafo, complex<double> c,
		complex<double> s, int i, int j);

	/**
	 * \brief Calculate Jacobi-Angles c, s for given elemtents i, j
	 * @param c Givens matrix sin element
	 * @param s Givens matrix sin element
	 * @param i Target index in A
	 * @param j Target index in A
	 * @param A Set of matrices
	 */
	void calculateAngles(complex<double>& c, complex<double>& s,
		int i, int j, const vector<Matrixcd>& A);

	/**
	 * \brief Build the Matrix G, that is required to Calculate Jacobi-angles
	 * @param i Target index in A
	 * @param j Target index in A
	 * @param A Matrices for which G is built
	 * @return G-matrix
	 */
	Matrixcd buildGMatrix(int i, int j, const vector<Matrixcd>& A);

	/**
	 * \brief Weight matrices, i.e. x_w = 0.5 (AX + XA)
	 * @param A Matrices to be weighted
	 * @param W Weighting matrix
	 */
	void weightMatrices(vector<Matrixcd>& A, const Matrixcd& W);

	/**
	 * \brief Build hessian from precalculated objects
	 */
	Matrixd rfoBuildHessian(const Matrixd& preHessian,
		const Vectord& grad, double a);

	/**
	 * \brief Interpret a Vectord(2) as a complex<double>
	 * @param vec Vectod(2)
	 * @return complex number
	 */
	complex<double> interpretComplex(const Vectord& vec);

	/**
	 * \brief Change of the Diagonals of a matrix under Givens-rotation
	 * @param A Matrices to be rotated
	 * @param p target index
	 * @param q target index
	 * @param c cos of alpha_ in givens matrix
	 * @param s sin of alpha_ in givens matrix
	 * @return Change of diagonality-measure
	 */
	pair<double,double> rotatedDiagonals(const Matrixcd& A, int p, int q,
		complex<double> c, complex<double> s);

	/**
	 * \brief rotate a Matrix
	 * @param A Matrix that is rotated
	 * @param p target index
	 * @param q target index
	 * @param c cos(alpha_) in Givens matrix
	 * @param s sin(alpha_) in Givens matrix
	 * @return Rotated matrix A
	 */
	Matrixcd rotate(const Matrixcd& A,
		int p, int q, complex<double> c, complex<double> s);
};

