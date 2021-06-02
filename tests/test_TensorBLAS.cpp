//
// Created by Roman Ellerbrock on 5/25/21.
//

#include "Core/TensorBLAS.h"
#include <UnitTest++/UnitTest++.h>
#include "Core/Tensor_Extension.h"


using namespace std;

SUITE (Tensor) {
	double eps = 1e-7;

	typedef complex<double> cd;

	class Factory {
	public:
		Factory()
			: gen_(129123), dim_(16) {
			S_ = Create(dim_);
			TensorShape shape3({dim_, dim_, dim_});
			A_ = Tensorcd(shape3);
			B_ = Tensorcd(shape3);
			A2_ = Tensorcd(shape3);
			Tensor_Extension::generate(A_, gen_);
			Tensor_Extension::generate(B_, gen_);
		}

		Matrixcd Create(size_t dim) {
			Matrixcd x(dim, dim);
			Tensor_Extension::generate(x, gen_);
			return x;
		}

		size_t dim_;
		Matrixcd S_;
		Tensorcd A_, A2_, B_, A2sol_;
		mt19937 gen_;
	};

	TEST_FIXTURE (Factory, Transpose) {
		Matrixcd ST1(dim_, dim_);
		Matrixcd ST2(dim_, dim_);
		transpose(&ST1[0], &S_[0], dim_, dim_);
		transpose2<cd, 4>(&ST2[0], &S_[0], dim_, dim_);
			CHECK_CLOSE(0., residual(ST1, ST2), eps);

		/// === non even dim ===
		size_t d = 23;
		auto X = Create(d);
		Matrixcd XT1(d, d);
		Matrixcd XT2(d, d);
		transpose(&XT1[0], &X[0], d, d);
		transpose2<cd, 4>(&XT2[0], &X[0], d, d);
			CHECK_CLOSE(0., residual(XT1, XT2), eps);
	}

	/// ===== matrixTensor ================================

	TEST_FIXTURE (Factory, matrixTensorT1) {
		/// Test matrixTensor routine for square matrices
		/// First, the correct solution
		A2sol_ = matrixTensor(S_, A_, 1);

		/// apply on mode 1
		Tensorcd workA(A2_.shape());
		matrixTensor2(A2_, S_, A_, workA, dim_, dim_, dim_, dim_, 1);
			CHECK_CLOSE(0., residual(A2_, A2sol_), eps);

		/// apply on mode 0
		A2sol_ = matrixTensor(S_, A_, 0);
		matrixTensor2(A2_, S_, A_, workA, 1, dim_, dim_, dim_ * dim_, 1);
			CHECK_CLOSE(0., residual(A2_, A2sol_), eps);

		/// apply on mode 2
		A2sol_ = matrixTensor(S_, A_, 2);
		matrixTensor2(A2_, S_, A_, workA, dim_ * dim_, dim_, dim_, 1, 1);
			CHECK_CLOSE(0., residual(A2_, A2sol_), eps);
	}

	TEST_FIXTURE (Factory, matrixTensorRectMode0) {
		Matrixcd Srec(2 * dim_, dim_);
		TensorShape shape({2 * dim_, dim_, dim_});
		Tensorcd A2(shape);
		Tensorcd A2sol(shape);
		Tensorcd workA = A2;
		matrixTensor2(A2, Srec, A_, workA, 1, 2 * dim_, dim_, dim_ * dim_, 1);
		matrixTensor1(A2sol, Srec, A_, 1, 2 * dim_, dim_, dim_ * dim_, 1);
			CHECK_CLOSE(0., residual(A2, A2sol), eps);
		matrixTensorBLAS(A2, Srec, A_, 0);
			CHECK_CLOSE(0., residual(A2, A2sol), eps);
	}

	TEST_FIXTURE (Factory, matrixTensorRectMode1) {
		Matrixcd Srec(2 * dim_, dim_);
		TensorShape shape({dim_, 2 * dim_, dim_});
		Tensorcd A2(shape);
		Tensorcd A2sol(shape);
		Tensorcd workA = A2;
		matrixTensor2(A2, Srec, A_, workA, dim_, 2 * dim_, dim_, dim_, 1);
		matrixTensor1(A2sol, Srec, A_, dim_, 2 * dim_, dim_, dim_, 1);
			CHECK_CLOSE(0., residual(A2, A2sol), eps);

		matrixTensorBLAS(A2, Srec, A_, 1);
			CHECK_CLOSE(0., residual(A2, A2sol), eps);
	}

	TEST_FIXTURE (Factory, matrixTensorRectMode2) {
		Matrixcd Srec(2 * dim_, dim_);
		TensorShape shape({dim_, dim_, 2 * dim_});
		Tensorcd A2(shape);
		Tensorcd A2sol(shape);
		Tensorcd workA = A2;
		matrixTensor2(A2, Srec, A_, workA, dim_ * dim_, 2 * dim_, dim_, 1, 1);
		matrixTensor1(A2sol, Srec, A_, dim_ * dim_, 2 * dim_, dim_, 1, 1);
			CHECK_CLOSE(0., residual(A2, A2sol), eps);
		matrixTensorBLAS(A2, Srec, A_, 2);
			CHECK_CLOSE(0., residual(A2, A2sol), eps);
	}

	/// ===== Contraction ================================

	TEST_FIXTURE (Factory, contractionT1) {
		Matrixcd h(dim_, dim_);

		auto hSol = contraction(A_, B_, 0);
		contractionBLAS(h, A_, B_, 0);
			CHECK_CLOSE(0., residual(hSol, h), eps);

		hSol = contraction(A_, B_, 1);
		contractionBLAS(h, A_, B_, 1);
			CHECK_CLOSE(0., residual(hSol, h), eps);

		hSol = contraction(A_, B_, 2);
		contractionBLAS(h, A_, B_, 2);
			CHECK_CLOSE(0., residual(hSol, h), eps);
	}

	TEST_FIXTURE (Factory, contractionT2mode0) {
		TensorShape shapeB({2 * dim_, dim_, dim_});
		Matrixcd h(dim_, 2 * dim_);
		Matrixcd hSol(dim_, 2 * dim_);
		Tensorcd Brec(shapeB);
		contraction1(hSol, A_, Brec, 1, dim_, 2 * dim_, dim_ * dim_, 0);
		Tensorcd workA(A_);
		Tensorcd workB(Brec);
		contraction2(h, A_, Brec, workA, workB, 1, dim_, 2 * dim_, dim_ * dim_, 0);
			CHECK_CLOSE(0., residual(h, hSol), eps);
	}

	TEST_FIXTURE (Factory, contractionT2mode1) {
		TensorShape shapeB({2 * dim_, dim_, dim_});
		Matrixcd h(dim_, 2 * dim_);
		Matrixcd hSol(dim_, 2 * dim_);
		Tensorcd Brec(shapeB);
		contraction1(hSol, A_, Brec, dim_, dim_, 2 * dim_, dim_, 0);
		Tensorcd workA(A_);
		Tensorcd workB(Brec);
		contraction2(h, A_, Brec, workA, workB, dim_, dim_, 2 * dim_, dim_, 0);
			CHECK_CLOSE(0., residual(h, hSol), eps);
	}

	TEST_FIXTURE (Factory, contractionT2mode2) {
		TensorShape shapeB({2 * dim_, dim_, dim_});
		Matrixcd h(dim_, 2 * dim_);
		Matrixcd hSol(dim_, 2 * dim_);
		Tensorcd Brec(shapeB);
		contraction1(hSol, A_, Brec, dim_, dim_ * dim_, 2 * dim_, 1, 0);
		Tensorcd workA(A_);
		Tensorcd workB(Brec);
		contraction2(h, A_, Brec, workA, workB, dim_ * dim_, dim_, 2 * dim_, 1, 0);
			CHECK_CLOSE(0., residual(h, hSol), eps);
	}

}

