#include "Util/SimultaneousDiagonalization.h"
#include "Tensor/TensorLapack.h"

namespace simultaneousDiagonalization {

	void calculate(vector<Matrixcd>& A,
		Matrixcd& trafo, double eps) {
		bool converged = false;
		int iter = 0;
		int maxiter = 100;

		// Initialize the Diagonalization by rotating to the diagonal
		// representation of one of the matrices in A. This
		// avoids stationary points during the optimization process.
		trafo.zero();
		for (size_t i = 0; i < trafo.shape_.front(); i++)
			trafo(i, i) = 1.;

		// This rotates to the eigenspace of the first matrix
//		initialTransformation(A, trafo);

		// Iterate Jacobirotations until a converged result is reached
		// Measure off-diagonal norm
		double delta = measureDiagonality(A);
		double delta_off = measureOffDiagonals(A);
		cout << "Start : " << delta << "\t" << delta_off << endl;
		while (!converged && iter < maxiter) {
			// Rotation circle over all elements
			A.front().print();
			cout << endl;
			jacobiRotations(A, trafo);
			A.front().print();

			// Measure off-diagonal norm
			delta = measureDiagonality(A);
			delta_off = measureOffDiagonals(A);

			// Check convergence
			iter++;
			cout << iter << " : " << delta << "\t" << delta_off << endl;

			if (delta < eps) { converged = true; }
		}
	}

	void jacobiRotations(vector<Matrixcd>& A,
		Matrixcd& trafo) {
		// Angles for Givens-Rotation
		complex<double> c, s = 0;

		// Swipe over the matrix-dimensions and perform jacobi-rotations
		size_t dim = trafo.shape_[0];
		for (size_t i = 0; i < dim - 1; i++) {
			for (size_t j = i + 1; j < dim; j++) {
				// Calculate Angles c and s for the elements i and j
				calculateAngles(c, s, i, j, A);

				assert(abs(1. - abs(c) * abs(c) - abs(s) * abs(s)) < 1E-10);

				// Perform the Givens-Rotation with angles c and s
				rotateMatrices(A, c, s, i, j);

				// Rotate the Transformation-Matrix
				givensTrafoRotation(trafo, c, s, i, j);
			}
		}
	}

	double measureOffDiagonals(const vector<Matrixcd>& A) {
		// Measure the norm of off-diagonal elements
		double eps = 0;

		for (size_t k = 0; k < A.size(); k++) {
			Matrixcd B(A[k]);

			for (int n = 0; n < B.shape_[0]; n++)
				B(n, n) = 0;

			B = B * B;
			eps += real(trace(B));
		}

		return sqrt(eps);
	}

	double measureDiagonality(vector<Matrixcd>& A) {
		// Measure the norm of off-diagonal elements
		double eps = 0;

		for (size_t k = 0; k < A.size(); k++) {
			Matrixcd& B = A[k];

			size_t dim = B.shape_.front();

			for (size_t n = 0; n < dim; n++)
				for (size_t m = 0; m < dim; m++)
					if (m != n) {
						eps += pow(abs(B(m, n)), 2);
					}
		}

		return sqrt(eps);
	}

	void initialTransformation(vector<Matrixcd>& A,
		Matrixcd& trafo) {
		Matrixcd& B = A[0];
//		Vectord ev(B.shape_.front());
		auto x = heev(B);
		trafo = x.U();
		const auto& ev = x.ev();
//		B.cDiag(trafo, ev);

		for (size_t k = 0; k < A.size(); k++) {
			A[k] = unitarySimilarityTrafo(A[k], trafo);
		}
	}
}
