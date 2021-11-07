//
// Created by Roman Ellerbrock on 11/7/21.
//

#ifndef TENSORFORTRAN_H
#define TENSORFORTRAN_H

extern "C" {
// subroutine matvec (mulpsi, psi, matrix, a, b, c, add)
// subroutine rhomat (bra,ket,matrix,a,b,c)
void matvec_(double *C, double *B, double *mat,
	int *a, int *b, int *c, int *add);
void ctmatvec_(double *C, double *B, double *mat,
	int *a, int *b, int *c, int *add);
void rmatvec_(double *C, double *B, double *mat,
	int *a, int *b, int *c, int *add);
void rhomat_(double *Bra, double *Ket, double *M,
	int *a, int *b, int *c);
}

template<typename T>
void contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind) {

	typedef complex<double> cd;
	typedef double d;
	if constexpr(is_same<T, cd>::value) {
		if (active1 == active2) {

			int a = active1;
			int b = before;
			int c = behind;

			rhomat_((double *) &A[0], (double *) &B[0], (double *) &S[0],
				&a, &b, &c);
			return;
		}
	}

	contraction1(S, A, B, before, active1, active2, behind, true);
}

template<typename T>
void contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind) {

	typedef complex<double> cd;
	typedef double d;
	if constexpr(is_same<T, cd>::value) {
		if (active1 == active2) {

			int a = active1;
			int b = before;
			int c = behind;

			rhomat_((double *) &A[0], (double *) &B[0], (double *) &S[0],
				&a, &b, &c);
			return;
		}
	}

	contraction1(S, A, B, before, active1, active2, behind, true);
}

#endif //TENSORFORTRAN_H
