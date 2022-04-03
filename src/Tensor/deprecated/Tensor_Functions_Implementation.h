#include "Tensor_Functions.h"

template<typename T>
Matrix<T> dotProduct(const Tensor<T>& A, const Tensor<T>& B) {
	size_t nmax = B.shape().lastDimension();
	size_t mmax = A.shape().lastDimension();
	Matrix<T> S(mmax, nmax);
	contraction(S, A, B, A.shape().lastIdx());
	return S;
}

template<typename T, typename U>
void contraction1(Matrix<U>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {

	if (zero) { h.zero(); }

	for (size_t a = 0; a < A; ++a) {
		for (size_t b = 0; b < B; ++b) {
			for (size_t b2 = 0; b2 < B2; ++b2) {
				for (size_t c = 0; c < C; ++c) {
					//					h(b, b2) += conj(bra(a, b, c)) * ket(a, b2, c);
					h(b, b2) += conj(bra[a + b * A + c * A * B]) * ket(a + b2 * A + c * A * B2);
				}
			}
		}
	}
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
Matrix<T> contraction(const Tensor<T>& A, const Tensor<T>& B, size_t k) {
	const TensorShape& tdim_a(A.shape());
	const TensorShape& tdim_b(B.shape());
	assert(k < tdim_a.order());
	assert(k < tdim_b.order());
	size_t active1 = tdim_a[k];
	size_t active2 = tdim_b[k];
	Matrix<T> S(active1, active2);
	contraction(S, A, B, k);
	return S;
}

template<typename T>
void contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B, size_t k, bool zero) {
	const TensorShape& tdim_a(A.shape());
	const TensorShape& tdim_b(B.shape());
	assert(k < tdim_a.order());
	assert(k < tdim_b.order());
	size_t before = tdim_a.before(k);
	size_t after = tdim_a.after(k);
	assert(tdim_b.before(k) == before);
	assert(tdim_b.after(k) == after);
	size_t active1 = tdim_a[k];
	size_t active2 = tdim_b[k];
	assert(tdim_a.totalDimension() / active1 == tdim_b.totalDimension() / active2);
	if (zero) { S.zero(); }
	contraction(S, A, B, before, active1, active2, after);
}

template<typename T, typename U>
void matrixTensor1(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero) {

	if (zero) { C.zero(); }

	size_t dimafter = active * before;
	size_t dimafterC = activeC * before;
	for (size_t aft = 0; aft < after; ++aft) {
		for (size_t act = 0; act < active; ++act) {
			for (size_t actC = 0; actC < activeC; ++actC) {
				for (size_t bef = 0; bef < before; ++bef) {
					C[dimafterC * aft + actC * before + bef] +=
						h[act * activeC + actC] * B[dimafter * aft + act * before + bef];
				}
			}
		}
	}
}

template<typename T, typename U>
void matrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero) {
	// Null the result tensor if flag is set to "true"
	int add = !zero;
	int a = activeB;
	int b = before;
	int c = after;
	typedef complex<double> cd;
	typedef double d;

	if (activeB == activeC) {
		if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
			matvec_((double *) &C[0], (double *) &B[0], (double *) &A[0],
				&a, &b, &c, &add);
			return;
		} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
			rmatvec_((double *) &C[0], (double *) &B[0], (double *) &A[0],
				&a, &b, &c, &add);
			return;
		}
	}

	matrixTensor1(C, A, B, before, activeC, activeB, after, zero);
}

template<typename T, typename U>
void tMatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero) {
	// Null the result tensor if flag is set to "true"

	int add = !zero;
	int a = activeB;
	int b = before;
	int c = after;
	typedef complex<double> cd;
	typedef double d;
	if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
		ctmatvec_((double *) &C[0], (double *) &B[0], (double *) &A[0],
			&a, &b, &c, &add);
		//	} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
		//		rmatvec_((double*)&C[0], (double*)&B[0], (double*)&A[0],
		//			&a, &b, &c, &add);
	} else {
		if (zero) { C.zero(); }
		size_t actbefB = activeB * before;
		size_t actbefC = activeC * before;
		size_t Cidx = 0;
		size_t Bidx = 0;
		size_t Aidx = 0;
		size_t kpreidxB = 0;
		size_t kpreidxC = 0;
		size_t Bpreidx = 0;
		size_t Cpreidx = 0;

		if (before == 1) {
			//#pragma omp parallel for
			for (size_t k = 0; k < after; ++k) {
				kpreidxB = k * actbefB;
				kpreidxC = k * actbefC;
				for (size_t l = 0; l < activeB; ++l) {
					Bidx = l + kpreidxB;
					for (size_t j = 0; j < activeC; ++j) {
						Cidx = j + kpreidxC;
						Aidx = j * activeB + l;
						C[Cidx] += conj(A[Aidx]) * B[Bidx];
						//					C[Cidx] += conj(A[l, j]) * B[Bidx];
					}
				}
			}
		} else {
			//#pragma omp parallel for
			for (size_t k = 0; k < after; ++k) {
				kpreidxB = k * actbefB;
				kpreidxC = k * actbefC;
				for (size_t l = 0; l < activeB; ++l) {
					Bpreidx = l * before + kpreidxB;
					for (size_t j = 0; j < activeC; ++j) {
						Aidx = j * activeB + l;
						Cpreidx = j * before + kpreidxC;
						for (size_t i = 0; i < before; ++i) {
							Cidx = Cpreidx + i;
							Bidx = Bpreidx + i;
							C[Cidx] += conj(A[Aidx]) * B[Bidx];
							//						C[Cidx] += conj(A[l, j]) * B[Bidx];
						}
					}
				}
			}
		}
	}
}

template<typename T, typename U>
void matrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {
	TensorShape tdim(B.shape());
	TensorShape tdimC(C.shape());

	size_t after = tdim.after(mode);
	size_t before = tdim.before(mode);
	size_t active1 = A.dim1();
	size_t active2 = A.dim2();

	assert(mode < tdim.order());
	assert(A.dim2() == tdim[mode]);
	assert(A.dim1() == tdimC[mode]);

	matrixTensor(C, A, B, before, active1, active2, after, zero);
}

template<typename T, typename U>
Tensor<T> matrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode) {
	const TensorShape& tdim(B.shape());
	assert(mode < tdim.order());

	if (A.dim1() == A.dim2()) {
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t active = tdim[mode];
		size_t before = tdim.before(mode);
		matrixTensor(C, A, B, before, active, active, after, false);
		return C;
	} else {
		TensorShape tdim(B.shape());
		size_t active1 = A.dim1();
		size_t active2 = A.dim2();
		tdim = replaceDimension(tdim, mode, active1);
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t before = tdim.before(mode);
		assert(active1 == C.shape()[mode]);
		assert(active2 == B.shape()[mode]);
		cout << "non-quadratic mattensor implemented but tested only once so far.\n";
		matrixTensor(C, A, B, before, active1, active2, after, false);
		cout << "done.\n";
		return C;
	}
}

template<typename T, typename U>
void tensorMatrix(Tensor<T>& C, const Tensor<T>& B, const Matrix<U>& A, size_t mode, bool zero) {
	tensorMatrix(C, B, A.transpose(), mode, zero);
}

template<typename T, typename U>
Tensor<T> tensorMatrix(const Tensor<T>& B, const Matrix<U>& A, size_t mode) {
	return matrixTensor(A.transpose(), B, mode);
}

template<typename T, typename U>
Tensor<T> tMatrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode) {
	/// @TODO: remove this function; replace by tensorMatrix
	const TensorShape& tdim(B.shape());
	assert(mode < tdim.order());
	assert(mode >= 0);
	assert(A.dim1() == B.shape()[mode]);

	if (A.dim1() == A.dim2()) {
		Tensor<T> C(tdim);
		size_t after = tdim.after(mode);
		size_t active = tdim[mode];
		size_t before = tdim.before(mode);
		tMatrixTensor(C, A, B, before, active, active, after, false);
		return C;
	} else {
		size_t activeC = A.dim2();
		size_t activeB = A.dim1();
		TensorShape tdim(B.shape());
		tdim = replaceDimension(tdim, mode, A.dim2());
		size_t after = tdim.after(mode);
		size_t before = tdim.before(mode);
		Tensor<T> C(tdim);
		cout << "non-quadratic mattensor implemented but not tested, yet.\n";
		tMatrixTensor(C, A, B, before, activeC, activeB, after, false);
		getchar();
		return C;
	}
}

template<typename T, typename U>
void multStateAB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, bool zero) {
	const TensorShape& tdimB(B.shape());
	const TensorShape& tdimC(C.shape());

	const size_t before = tdimB.lastBefore();
	const size_t active1 = tdimB.lastDimension();
	const size_t active2 = tdimC.lastDimension();
	const size_t after = 1;

	assert(A.dim2() == active1);
	assert(A.dim1() == active2);
	assert(before == tdimC.lastBefore());

	matrixTensor(C, A, B, before, active1, active2, after, zero);
}

template<typename T, typename U>
Tensor<T> multStateAB(const Matrix<U>& A, const Tensor<T>& B) {
	const TensorShape& tdim_b(B.shape());
	size_t ntensor = tdim_b.lastDimension();
	assert(A.dim2() == ntensor);

	TensorShape tdim_c(tdim_b);
	tdim_c.setDimension(A.dim1(), tdim_c.lastIdx());
	Tensor<T> C(tdim_c);
	multStateAB(C, A, B);
	return C;
}

template<typename T, typename U>
void multStateArTB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B) {
	const TensorShape& tdim(B.shape());
	size_t dimpart = tdim.lastBefore();
	size_t ntensor = tdim.lastDimension();
	for (size_t n = 0; n < ntensor; n++) {
		size_t B_idx = n * dimpart;
		for (size_t m = 0; m < ntensor; m++) {
			size_t C_idx = m * dimpart;
			size_t A_idx = m * ntensor;
			for (size_t i = 0; i < dimpart; i++) {
				/// C(i, m) += A(n, m) * B(i, n);
				C[C_idx + i] += A[A_idx + n] * B[B_idx + i];
			}
		}
	}
}

template<typename T, typename U>
Tensor<T> multStateArTB(const Matrix<U>& A, const Tensor<T>& B) {
	const TensorShape& tdim(B.shape());
	size_t dimpart = tdim.lastBefore();
	size_t ntensor = tdim.lastDimension();
	assert(A.dim1() == A.dim2());
	assert(A.dim2() == ntensor);

	Tensor<T> C(tdim);
	multStateArTB(C, A, B);

	return C;
}

template<typename T>
double residual(Tensor<T> A, const Tensor<T>& B) {
	A -= B;
	auto S = dotProduct(A, A);
	return S.frobeniusNorm();
}


template<typename T>
Matrix<T> toMatrix(const Tensor<T>& A) {
	const TensorShape& shape = A.shape();
	size_t diml = shape.lastBefore();
	size_t dimr = shape.lastDimension();
	Matrix<T> B(diml, dimr);
	for (size_t j = 0; j < dimr; ++j) {
		for (size_t i = 0; i < diml; ++i) {
			B(i, j) = A(i, j);
		}
	}
	return B;
}

template<typename T>
Matrix<T> moveToMatrix(Tensor<T>& A) {
	const TensorShape& shape = A.shape();
	size_t diml = shape.lastBefore();
	size_t dimr = shape.lastDimension();
	A.ownership_ = false;
	Matrix<T> B(diml, dimr, A.coeffs_, true, false);
	return B;
}

template<typename T>
Matrix<T> toMatrix(const Tensor<T>& A, size_t mode) {
	const TensorShape& shape = A.shape();
	size_t dimbef = shape.before(mode);
	size_t dimaft = shape.after(mode);
	size_t diml = dimbef * dimaft;
	size_t dimr = shape[mode];
	Matrix<T> B(diml, dimr);
	TensorShape tmp({dimbef, dimr, dimaft});
	for (size_t bef = 0; bef < dimbef; ++bef) {
		for (size_t a = 0; a < dimr; ++a) {
			for (size_t aft = 0; aft < dimaft; ++aft) {
				size_t idxl = bef + aft * dimbef;
				size_t I = indexMapping({bef, a, aft}, tmp);
				B(idxl, a) = A(I);
			}
		}
	}
	return B;
}

template<typename T>
Tensor<T> toTensor(const Matrix<T>& B, const TensorShape& shape, size_t mode) {
	Tensor<T> A(shape);

	size_t dimbef = shape.before(mode);
	size_t dimaft = shape.after(mode);
	size_t dimr = shape[mode];
	TensorShape tmp({dimbef, dimr, dimaft});
	for (size_t bef = 0; bef < dimbef; ++bef) {
		for (size_t a = 0; a < dimr; ++a) {
			for (size_t aft = 0; aft < dimaft; ++aft) {
				size_t idxl = bef + aft * dimbef;
				size_t I = indexMapping({bef, a, aft}, tmp);
				A(I) = B(idxl, a);
			}
		}
	}
	return A;
}

template<typename T>
Tensor<T> toTensor(const Matrix<T>& B) {
	TensorShape shape({B.dim1(), B.dim2()});
	Tensor<T> A(shape);
	for (size_t j = 0; j < B.dim2(); ++j) {
		for (size_t i = 0; i < B.dim1(); ++i) {
			A(i, j) = B(i, j);
		}
	}
	return A;
}

template<typename T>
Tensor<T> moveToTensor(Matrix<T>& B) {
	TensorShape shape({B.dim1(), B.dim2()});
	Tensor<T> A(shape, &B[0], true, false);
	return A;
}




//Projects B on A
template<typename T>
Tensor<T> project(const Tensor<T>& A,
	const Tensor<T>& B) {
	//calculates the overlap of A with it self
	Tensor<T> Aperp(A);
	//	GramSchmidt(Aperp);
	const Matrix<T> overlap = Aperp.dotProduct(Aperp);

	//invert the overlap
	const Matrix<T> inverse_operlap = overlap.cInv();

	//calculate the scalar product of A and B
	const Matrix<T> dotproduct = Aperp.dotProduct(B);

	//multiply the scalar product and the inverse_operlap
	const Matrix<T> product = inverse_operlap * dotproduct;

	return multStateArTB(product, Aperp);
	//	return multStateArTB(dotproduct, Aperp);
}

/*! \brief Project B out of A, i.e. Anew = (1-P_B)*A
 *
 * This routine takes a Tensor A and makes it orthogonal to B.
 * It can be written as A_new = (1 - P_B) A, where P_B is the
 * projector onto B.
 */
template<typename T>
Tensor<T> projectOut(const Tensor<T>& A,
	const Tensor<T>& B) {
	Tensor<T> projector = project(B, A);
	Tensor<T> perp_A(A);
	const TensorShape& tdim = A.shape();
	for (size_t i = 0; i < tdim.totalDimension(); ++i) {
		perp_A(i) -= projector(i);
	}
	return perp_A;
}

//Projects B on A
template<typename T>
Tensor<complex<double> > projectOrthogonal(const Tensor<complex<double> >& A,
	const Tensor<T>& B) {
	// calculate the scalar product of A and B
	const Matrix<complex<double> > dotproduct = dotProduct(A, B);

	return multStateArTB(dotproduct, A);
}





tuple<Tensorcd, Matrixcd, Vectord> SVD(const Tensorcd& A) {
	const TensorShape& tdim = A.shape();
	size_t dimpart = tdim.lastBefore();
	size_t ntensor = tdim.lastDimension();

	using namespace Eigen;
	MatrixXcd Am = Eigen::Map<MatrixXcd>((complex<double> *) &A(0), dimpart, ntensor);
	JacobiSVD<MatrixXcd> svd(Am, ComputeThinU | ComputeThinV);

	Tensorcd U(tdim);
	auto u_mat = svd.matrixU();
	for (size_t i = 0; i < dimpart; ++i) {
		for (size_t n = 0; n < ntensor; ++n) {
			U(i, n) = u_mat(i, n);
		}
	}

	auto v_mat = svd.matrixV();
	Matrixcd V(ntensor, ntensor);
	for (size_t i = 0; i < ntensor; ++i) {
		for (size_t n = 0; n < ntensor; ++n) {
			V(i, n) = v_mat(i, n);
		}
	}

	auto sigma_e = svd.singularValues();
	Vectord sigma(ntensor);
	for (size_t i = 0; i < ntensor; ++i) {
		sigma(i) = sigma_e(i);
	}

	return tuple<Tensorcd, Matrixcd, Vectord>(U, V, sigma);
}

tuple<Matrixcd, Matrixcd, Vectord> SVD(const Matrixcd& A) {
	size_t dim1 = A.dim1();
	size_t dim2 = A.dim2();

	using namespace Eigen;
	MatrixXcd Am = Eigen::Map<MatrixXcd>((complex<double> *) &A(0, 0), dim1, dim2);
	JacobiSVD<MatrixXcd> svd(Am, ComputeThinU | ComputeThinV);

	auto u_mat = svd.matrixU();
	dim1 = u_mat.rows();
	dim2 = u_mat.cols();
	Matrixcd U(dim1, dim2);
	for (size_t i = 0; i < dim1; ++i) {
		for (size_t n = 0; n < dim2; ++n) {
			U(i, n) = u_mat(i, n);
		}
	}

	auto v_mat = svd.matrixV();
	dim1 = v_mat.rows();
	dim2 = v_mat.cols();
	Matrixcd V(dim1, dim2);
	for (size_t i = 0; i < dim1; ++i) {
		for (size_t n = 0; n < dim2; ++n) {
			V(i, n) = v_mat(i, n);
		}
	}

	auto sigma_e = svd.singularValues();
	size_t mdim = min(dim1, dim2);
	Vectord sigma(mdim);
	for (size_t i = 0; i < mdim; ++i) {
		sigma(i) = sigma_e(i);
	}

	return tuple<Matrixcd, Matrixcd, Vectord>(U, V, sigma);
}

Tensorcd normalize(Tensorcd A, size_t k, double eps) {
	/**
	 * \brief Normalize a tensor for the k-th leg while preserving its representation
	 */
	auto a = toMatrix(A, k);
	auto x = svd(a);
	Matrixcd& U = get<0>(x);
	const Vectord& sigma = get<2>(x);

	/// fill
	mt19937 gen(time(nullptr));
	TensorShape Ushape({U.dim1(), U.dim2()});
	Tensorcd rand(Ushape);
	Tensor_Extension::generate(rand, gen);

	for (size_t j = 0; j < U.dim2(); ++j) {
		if (sigma(j) < eps) {
			for (size_t i = 0; i < U.dim1(); ++i) {
				U(i, j) = rand(i, j);
			}
		}
	}
	auto utens = toTensor(U);
	U = toMatrix(utens);
	gramSchmidt(utens);
	U = toMatrix(utens);

	const Matrixcd& V = get<1>(x);
	a = U * V.adjoint();
	A = toTensor(a, A.shape(), k);
	return A;
}

template<typename T>
Matrix<T> map(const Tensor<T>& A) {
	const TensorShape& tdim = A.shape();
	size_t ntensor = tdim.lastDimension();
	size_t dimpart = tdim.lastBefore();
	Matrix<T> M(dimpart, ntensor);
	for (size_t n = 0; n < ntensor; ++n) {
		for (size_t i = 0; i < dimpart; ++i) {
			M(i, n) = A(i, n);
		}
	}
	return M;
}
// @TODO:: Add mapping to and from Eigen

//////////////////////////////////////////////////////////////////////
/// Direct Sum + Product
//////////////////////////////////////////////////////////////////////

/*	template<typename T>
	Tensor<T> Merge(Tensor<T> A, const Tensor<T>& B) {
		// Merge two Tensors into one.
		const TensorShape& tdim1 = A.shape();
		const TensorShape& tdim2 = B.shape();
		size_t ntens1 = tdim1.lastDimension();
		size_t ntens2 = tdim2.lastDimension();
		A = A.AdjustStateDim(ntens1 + ntens2);
		for (size_t n = 0; n < ntens2; ++n) {
			for (size_t i = 0; i < tdim1.lastBefore(); ++i) {
				A(i, ntens1 + n) = B(i, n);
			}
		}
		return A;
	}*/

void shiftIndices(vector<size_t>& Ibreak, const TensorShape& shift,
	bool beforeLast, bool last) {
	if (beforeLast) {
		for (size_t k = 0; k < shift.lastIdx(); ++k) {
			Ibreak[k] += shift[k];
		}
	}
	if (last) {
		size_t idx = shift.lastIdx();
		Ibreak[idx] += shift[idx];
	}
}

TensorShape directSum(const TensorShape& A, const TensorShape& B,
	bool before, bool last) {
	assert(A.order() == B.order());
	vector<size_t> dims = A.dimensions();
	shiftIndices(dims, B, before, last);
	return TensorShape(dims);
}

template<typename T>
Tensor<T> directSum(const Tensor<T>& A, const Tensor<T>& B,
	bool before, bool last) {

	TensorShape Cshape = directSum(A.shape(), B.shape(), before, last);
	const TensorShape& Ashape = A.shape();
	const TensorShape& Bshape = B.shape();
	Tensor<T> C(Cshape);
	/// Place elements of A into C
	for (size_t I = 0; I < Ashape.totalDimension(); ++I) {
		auto Ibreak = indexMapping(I, Ashape);
		auto L = indexMapping(Ibreak, Cshape);
		C(L) = A(I);
	}
	/// Place elements of B into C
	for (size_t I = 0; I < Bshape.totalDimension(); ++I) {
		auto Ibreak = indexMapping(I, Bshape);
		shiftIndices(Ibreak, Ashape, before, last);
		auto L = indexMapping(Ibreak, Cshape);
		C(L) = B(I);
	}
	return C;
}

TensorShape DirectProduct(const TensorShape& A, const TensorShape& B) {
	assert(A.order() == B.order());
	auto dims = A.dimensions();
	for (size_t k = 0; k < A.order(); ++k) {
		dims[k] *= B[k];
	}
	return TensorShape(dims);
}

size_t mergeIndex(size_t I, size_t J, const TensorShape& A,
	const TensorShape& B, const TensorShape& C) {
	auto Ibreak = indexMapping(I, A);
	auto Jbreak = indexMapping(J, B);
	auto Lbreak(Ibreak);
	for (size_t k = 0; k < Ibreak.size(); ++k) {
		Lbreak[k] = Jbreak[k] * A[k] + Ibreak[k];
	}
	return indexMapping(Lbreak, C);
}

template<typename T>
Tensor<T> directProduct(const Tensor<T>& A, const Tensor<T>& B) {
	TensorShape shape = directProduct(A.shape(), B.shape());
	Tensor<T> C(shape);
	for (size_t J = 0; J < B.shape().totalDimension(); ++J) {
		for (size_t I = 0; I < A.shape().totalDimension(); ++I) {
			J = mergeIndex(I, J, A.shape(), B.shape());
			C(J) = A(I) * B(J);
		}
	}
	return C;
}

template<typename T>
Tensor<T> doubleHoleContraction(const Tensor<T>& A, const Tensor<T>& B,
	size_t k1, size_t k2) {

	const TensorShape& shape = A.shape();
	assert(k1 < shape.order());
	assert(k2 < shape.order());
	size_t dim1 = shape[k1];
	size_t dim2 = shape[k2];

	TensorShape nshape({dim1, dim2, dim1, dim2});
	Tensor<T> C(nshape);
	vector<size_t> Lbreak({0, 0, 0, 0});
	for (size_t I = 0; I < shape.totalDimension(); ++I) {
		auto Ibreak = indexMapping(I, shape);
		auto I2break = Ibreak;
		for (size_t l1 = 0; l1 < dim1; ++l1) {
			for (size_t l2 = 0; l2 < dim2; ++l2) {
				Lbreak[0] = Ibreak[k1];
				Lbreak[1] = Ibreak[k2];
				Lbreak[2] = l1;
				Lbreak[3] = l2;
				I2break[k1] = l1;
				I2break[k2] = l2;
				C(Lbreak) += conj(A(Ibreak)) * B(I2break);
			}
		}
	}
	return C;
}

//////////////////////////////////////////////////////////////////////
/// Random number routines for tensors
//////////////////////////////////////////////////////////////////////

/// Randomly occupy Tensors and Matrices
template<typename T>
void generateNormal(T *A, size_t n, mt19937& gen) {
	uniform_real_distribution<double> dist(-1., 1.);
	for (size_t i = 0; i < n; ++i) {
		A[i] = dist(gen);
	}
}

template<typename T>
void generate(Tensor<T>& A, mt19937& gen) {
	generateNormal(&A[0], A.shape().totalDimension(), gen);
}

template<typename T>
void generate(Matrix<T>& A, mt19937& gen) {
	generateNormal(&A[0], A.dim1() * A.dim2(), gen);
}

template<typename T>
void generate(Vector<T>& A, mt19937& gen) {
	generateNormal(&A[0], A.dim(), gen);
}

/* //////////////////////////////////////////////
 * Extension of the Tensor class
 *
 * The following functions are excluded from the
 * Tensor class to keep it slim.
 */

template<typename T>
void OuterProductAdd(Matrix<T>& M,
	const Tensor<T>& A, const Tensor<T>& B) {
	const TensorShape& tdim = A.shape();
	size_t dimpart = tdim.lastBefore();
	size_t ntensor = tdim.lastDimension();

#pragma omp parallel for
	for (size_t i = 0; i < dimpart; i++) {
		for (size_t j = 0; j < dimpart; j++) {
			for (size_t n = 0; n < ntensor; n++) {
				M(i, j) += A(i, n) * conj(B(j, n));
			}
		}
	}
}

template<typename T>
Matrix<T> outerProduct(const Tensor<T>& A, const Tensor<T>& B) {

	const TensorShape& tdim = A.shape();
	size_t dimpart = tdim.lastBefore();
	Matrix<T> M(dimpart, dimpart);
	OuterProductAdd(M, A, B);
	return M;
}

template<typename T>
void weightedOuterProductAdd(Matrix<T>& M, const Tensor<T>& A,
	const Tensor<T>& B, const Matrix<T>& rho) {
	Tensor<T> mA = multStateArTB(rho, A);
	OuterProductAdd(M, mA, B);
}

template<typename T>
Matrix<T> weightedOuterProduct(const Tensor<T>& A, const Tensor<T>& B,
	const Matrix<T>& m) {
	Tensor<T> mA = multStateAB(m, A);
//		Tensor<T> mB = multStateAB(m, B);
	return outerProduct(mA, B);
}

