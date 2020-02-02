#pragma once
#include "FactorMatrix.h"
#include "Vector.h"
#include "JacobiRotationFramework.h"

using namespace JacobiRotationFramework;

namespace WeightedSimultaneousDiagonalization
/**
 * \namespace WeightedSimultaneousDiagonalization
 * \ingroup Core
 * \brief Perform matrix-weighted Joint Diagonalization
 */
{
	/** \brief  Measure of the diagonality in WSD
	 *
	 * @param A Matrices to be checked
	 * @param W Weight matrix
	 * @return Measure of diagonality
	 */
	double MeasureWeightedDiagonality(const vector<FactorMatrixcd>& A,
			const FactorMatrixcd& W);

	/** \brief Calculate the weighted Simultaneous Diagonalization for the Matrices A with the weight W
	 *
	 * This performs a weighted simultaneous diagonalization of a set of matrices.
	 *
	 * @param Xs Matrices to diagonalize
	 * @param XXs Second moments of matrices (not used, can be empty)
	 * @param W Weight matrix
	 * @param trafo Output transformation
	 * @param eps Target accuracy
	 */
	void Calculate(vector<FactorMatrixcd>& Xs, vector<FactorMatrixcd> XXs, FactorMatrixcd& W,
		FactorMatrixcd& trafo, double eps);

	/// Quasi-Protected functions

	/**
	 * \brief Measure the weighted off-diagonality of a set of matrices
	 * @param Xws Weighted Matrices, i.e. 0.5(WX+XW)
	 * @param Xs Un-weighted matrices
	 * @param W Weight matrix
	 * @param trafo Transformation matrix
	 */
	double MeasureWeightedOffDiagonality(const vector<FactorMatrixcd>& Xws,
			const vector<FactorMatrixcd>& Xs, const FactorMatrixcd& W, const FactorMatrixcd& trafo);

	/**
	 * \brief Sweep over the whole matrix to optimize localization measure
	 * @param Xs
	 * @param XXs
	 * @param W
	 * @param trafo
	 */
	void WeightedJacobiRotations(vector<FactorMatrixcd>& Xs,
			vector<FactorMatrixcd>& XXs, FactorMatrixcd& W, FactorMatrixcd& trafo);

	/**
	 * \brief Calculate the optimal angles for WSD using Rational function optimizer
	 * @param c cos(alpha_); Element in Jacobi-matrix
	 * @param s sin(alpha_); Element in Jacobi-matrix
	 * @param i target index in X
	 * @param j target index in X
	 * @param X Matrices that should be diagonalized
	 * @param XXs Squared of X-matrices
	 * @param W Weight matrix
	 * @return error indicator
	 */
	int CalculateWeightedAngles(complex<double>& c,
			complex<double>& s, size_t i, size_t j,
			const vector<FactorMatrixcd>& X, const vector<FactorMatrixcd>& XXs,
			const FactorMatrixcd& W);

	/**
	 * \brief Calculate the change of the weighted-localization measure under givens rotation
	 * @param Xs Matrices to diagonalize
	 * @param XXs Squares of matrices that should be diagonalized
	 * @param W Weight matrices
	 * @param p target index in matrix
	 * @param q target index in matrix
	 * @param c element in Jacobi matrix
	 * @param s element in Jacobi matrix
	 * @return change of diagonality measure
	 */
	double WeightedJacobiLoc(const vector<FactorMatrixcd>& Xs,
			const vector<FactorMatrixcd>& XXs, const FactorMatrixcd& W,
			size_t p, size_t q, complex<double> c, complex<double> s);

	/**
	 * \brief Calculate derivatives (gradient and hessian) of the localization measurement
	 * @param grad Gradient of loc. measure
	 * @param preHessian hessian
	 * @param s_in complex-angle for rotation
	 * @param Xs matrices to diagonalize
	 * @param XXs squares of matrices to diagonalize
	 * @param W Weight matrix
	 * @param p target index in matrix
	 * @param q target index in matrix
	 * @param delta step length for numerical gradient
	 */
	void WeightedJacobiDerivatives(Vectord& grad, Matrixd& preHessian,
			complex<double> s_in, const vector<FactorMatrixcd>& Xs,
			const vector<FactorMatrixcd>& XXs, const FactorMatrixcd& W,
			size_t p, size_t q, double delta);
}

namespace WSD = WeightedSimultaneousDiagonalization;
