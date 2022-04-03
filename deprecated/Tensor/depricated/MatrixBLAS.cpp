//
// Created by Roman Ellerbrock on 6/4/21.
//

#include "Core/MatrixBLAS.h"
#include <cblas.h>
#include <lapacke.h>

void qrBLAS(Matrixcd& Q, const Matrixcd& A) {
	Matrixcd work(Q.dim1(), 1);
	Matrixcd tau(Q.dim1(), 1);
	Q = A;
	size_t m = A.dim1();
	size_t n = A.dim2();
	lapack_complex_double *a = reinterpret_cast<double _Complex *>(&Q[0]);
	lapack_complex_double *ta = reinterpret_cast<double _Complex *>(&tau[0]);
	lapack_int lda = m;
	lapack_int k = n;
	LAPACKE_zgeqrfp(LAPACK_COL_MAJOR, m, n, a, lda, ta);
	LAPACKE_zungqr(LAPACK_COL_MAJOR, m, n, k, a, lda, ta);
}

