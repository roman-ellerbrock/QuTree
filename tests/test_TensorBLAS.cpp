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


    TEST(onlyLegalContractions){

        // prepare tensors to be contracted
        const Tensorcd A{1,2,3,4};
        const Tensorcd B{1,4,5,2};
        const Tensorcd C{1,4,3,2,5};
        const Tensorcd D{1,2};

        // prepare index contracations
        std::vector<size_t> A_indices{0,1,3}; // legal
        std::vector<size_t> B_indices{0,3,1}; // also legal
        std::vector<size_t> C_indices{0,3,3}; // illegal
        std::vector<size_t> D_indices{0,1,3,4}; // illegal
        std::vector<size_t> E_indices{0,1}; // illegal
        std::vector<size_t> F_indices{2}; // illegal

        CHECK_EQUAL(true,is_contraction_legal(A,B,A_indices,B_indices));
        CHECK_EQUAL(true,is_contraction_legal(A,C,A_indices,B_indices));
        CHECK_EQUAL(true,is_contraction_legal(C,A,A_indices,B_indices));

        CHECK_EQUAL(false,is_contraction_legal(A,B,A_indices,C_indices));
        CHECK_EQUAL(false,is_contraction_legal(A,B,A_indices,E_indices));
        CHECK_EQUAL(false,is_contraction_legal(A,D,A_indices,C_indices));
        CHECK_EQUAL(false,is_contraction_legal(A,B,F_indices,F_indices));
    }

    TEST(generalTransposeOfFiveIndices){
        // prepare tensor
        const Tensor<double> A({2,3,4,5});
        for(int i = 0; i < A.shape().totalDimension(); ++i){
            A.coeffs_[i] = i;
        }

        Tensor<double> test({2,3,4,5});

        // check some transpose
        general_transpose_bd(test.coeffs_,A.coeffs_,2,3,4,5,1);

        // before: {1,0,2,3} is at position 85, now {1,3,2,0} is at position 27
        CHECK_EQUAL(85,A.coeffs_[85]);
        CHECK_EQUAL(85,test.coeffs_[27]);
    }

    TEST(generalTransposeOfGeneralTensor){
        // prepare tensor
        const Tensor<double> A({2,3,4,5,6,7,8});
        Tensor<double> res{A};

        general_transpose(res,A,6,3);

        CHECK_EQUAL(2, res.shape().dimensions()[0]);
        CHECK_EQUAL(3, res.shape().dimensions()[1]);
        CHECK_EQUAL(4, res.shape().dimensions()[2]);
        CHECK_EQUAL(8, res.shape().dimensions()[3]);
        CHECK_EQUAL(6, res.shape().dimensions()[4]);
        CHECK_EQUAL(7, res.shape().dimensions()[5]);
        CHECK_EQUAL(5, res.shape().dimensions()[6]);
    }

    TEST(generalTrasposeAndCorrectReorderingIsWorking){
        // prepare tensor
        const Tensor<double> A({2,3,4,5,6,7,8});
        Tensor<double> res{A};

        std::vector<size_t> new_form{6,0,3,2};

        general_transpose_to_order(res,A,new_form);

        CHECK_EQUAL(8, res.shape().dimensions()[0]);
        CHECK_EQUAL(2, res.shape().dimensions()[1]);
        CHECK_EQUAL(5, res.shape().dimensions()[2]);
        CHECK_EQUAL(4, res.shape().dimensions()[3]);
        CHECK_EQUAL(3, res.shape().dimensions()[4]);
        CHECK_EQUAL(6, res.shape().dimensions()[5]);
        CHECK_EQUAL(7, res.shape().dimensions()[6]);
    }

    TEST(MatrixMultiplicationAsGeneralTensorContraction){
        // first test: general matrix-matrix-multiplication
        Tensor<double> A{2,3};
        Tensor<double> B{3,4};
        Tensor<double> res{0};
        vector<size_t> contractionA{1};
        vector<size_t> contractionB{0};
        for(int i = 0; i < A.shape().totalDimension(); ++i){
            A.coeffs_[i] = i;
        }
        for(int i = 0; i < B.shape().totalDimension(); ++i){
            B.coeffs_[i] = i;
        }
        general_contraction(A,B,res,contractionA,contractionB);

        // test the result shape
        CHECK_EQUAL(2,res.shape().order());
        CHECK_EQUAL(2,res.shape().dimensions()[0]);
        CHECK_EQUAL(4,res.shape().dimensions()[1]);

        // test some of the result values
        CHECK_CLOSE(10.,res.coeffs_[0],1e-12); // (0,0)
        CHECK_CLOSE(13.,res.coeffs_[1],1e-12); // (1,0)
        CHECK_CLOSE(28.,res.coeffs_[2],1e-12); // (0,1)
        CHECK_CLOSE(40.,res.coeffs_[3],1e-12); // (1,1)
        CHECK_CLOSE(46.,res.coeffs_[4],1e-12); // (0,2)
        CHECK_CLOSE(67.,res.coeffs_[5],1e-12); // (1,2)
        CHECK_CLOSE(64.,res.coeffs_[6],1e-12); // (0,3)
        CHECK_CLOSE(94.,res.coeffs_[7],1e-12); // (1,3)

    }

    TEST(generalComplexTensorContraction){
        const Tensorcd A{1,2,3,4};
        const Tensorcd B{1,5,4,2};
        Tensorcd C{1,2,5,4};
        std::vector<size_t> A_indices{0,1,3};
        std::vector<size_t> B_indices{0,3,2};

        // prepare some demo data
        for(int i = 0; i < A.shape().totalDimension(); ++i){
            A.coeffs_[i] = std::complex<double>(i-1.,i+1.);
        }

        for(int i = 0; i < B.shape().totalDimension(); ++i){
            B.coeffs_[i] = std::complex<double>(i-1.,i+1.);
        }

        general_contraction(A,B,C,A_indices,B_indices);

        // result should be 2x3
        CHECK_EQUAL(3,C.shape().dimensions()[0]);
        CHECK_EQUAL(5,C.shape().dimensions()[1]);

        // check result values
                CHECK_CLOSE(   3356.0000000000000      ,C.coeffs_[           0 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[           0 ].imag(),1e-12);
                CHECK_CLOSE(   3916.0000000000000      ,C.coeffs_[           1 ].real(),1e-12);
                CHECK_CLOSE(  -96.000000000000000      ,C.coeffs_[           1 ].imag(),1e-12);
                CHECK_CLOSE(   4476.0000000000000      ,C.coeffs_[           2 ].real(),1e-12);
                CHECK_CLOSE(  -64.000000000000000      ,C.coeffs_[           2 ].imag(),1e-12);
                CHECK_CLOSE(   3508.0000000000000      ,C.coeffs_[           3 ].real(),1e-12);
                CHECK_CLOSE(  -144.00000000000000      ,C.coeffs_[           3 ].imag(),1e-12);
                CHECK_CLOSE(   4100.0000000000000      ,C.coeffs_[           4 ].real(),1e-12);
                CHECK_CLOSE(  -112.00000000000000      ,C.coeffs_[           4 ].imag(),1e-12);
                CHECK_CLOSE(   4692.0000000000000      ,C.coeffs_[           5 ].real(),1e-12);
                CHECK_CLOSE(  -80.000000000000000      ,C.coeffs_[           5 ].imag(),1e-12);
                CHECK_CLOSE(   3660.0000000000000      ,C.coeffs_[           6 ].real(),1e-12);
                CHECK_CLOSE(  -160.00000000000000      ,C.coeffs_[           6 ].imag(),1e-12);
                CHECK_CLOSE(   4284.0000000000000      ,C.coeffs_[           7 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[           7 ].imag(),1e-12);
                CHECK_CLOSE(   4908.0000000000000      ,C.coeffs_[           8 ].real(),1e-12);
                CHECK_CLOSE(  -96.000000000000000      ,C.coeffs_[           8 ].imag(),1e-12);
                CHECK_CLOSE(   3812.0000000000000      ,C.coeffs_[           9 ].real(),1e-12);
                CHECK_CLOSE(  -176.00000000000000      ,C.coeffs_[           9 ].imag(),1e-12);
                CHECK_CLOSE(   4468.0000000000000      ,C.coeffs_[          10 ].real(),1e-12);
                CHECK_CLOSE(  -144.00000000000000      ,C.coeffs_[          10 ].imag(),1e-12);
                CHECK_CLOSE(   5124.0000000000000      ,C.coeffs_[          11 ].real(),1e-12);
                CHECK_CLOSE(  -112.00000000000000      ,C.coeffs_[          11 ].imag(),1e-12);
                CHECK_CLOSE(   3964.0000000000000      ,C.coeffs_[          12 ].real(),1e-12);
                CHECK_CLOSE(  -192.00000000000000      ,C.coeffs_[          12 ].imag(),1e-12);
                CHECK_CLOSE(   4652.0000000000000      ,C.coeffs_[          13 ].real(),1e-12);
                CHECK_CLOSE(  -160.00000000000000      ,C.coeffs_[          13 ].imag(),1e-12);
                CHECK_CLOSE(   5340.0000000000000      ,C.coeffs_[          14 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[          14 ].imag(),1e-12);


    }

    TEST(generalComplexTensorContractionViaPairsWorks){
        const Tensorcd A{1,2,3,4};
        const Tensorcd B{1,5,4,2};
        Tensorcd C{1,2,5,4};

        const std::vector<std::pair<int,int>> contraction_pairs{{0,0},{1,3},{3,2}};

        // prepare some demo data
        for(int i = 0; i < A.shape().totalDimension(); ++i){
            A.coeffs_[i] = std::complex<double>(i-1.,i+1.);
        }

        for(int i = 0; i < B.shape().totalDimension(); ++i){
            B.coeffs_[i] = std::complex<double>(i-1.,i+1.);
        }

        general_contraction(A,B,C,contraction_pairs);

        // result should be 3x5
        CHECK_EQUAL(3,C.shape().dimensions()[0]);
        CHECK_EQUAL(5,C.shape().dimensions()[1]);

        // check result values
                CHECK_CLOSE(   3356.0000000000000      ,C.coeffs_[           0 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[           0 ].imag(),1e-12);
                CHECK_CLOSE(   3916.0000000000000      ,C.coeffs_[           1 ].real(),1e-12);
                CHECK_CLOSE(  -96.000000000000000      ,C.coeffs_[           1 ].imag(),1e-12);
                CHECK_CLOSE(   4476.0000000000000      ,C.coeffs_[           2 ].real(),1e-12);
                CHECK_CLOSE(  -64.000000000000000      ,C.coeffs_[           2 ].imag(),1e-12);
                CHECK_CLOSE(   3508.0000000000000      ,C.coeffs_[           3 ].real(),1e-12);
                CHECK_CLOSE(  -144.00000000000000      ,C.coeffs_[           3 ].imag(),1e-12);
                CHECK_CLOSE(   4100.0000000000000      ,C.coeffs_[           4 ].real(),1e-12);
                CHECK_CLOSE(  -112.00000000000000      ,C.coeffs_[           4 ].imag(),1e-12);
                CHECK_CLOSE(   4692.0000000000000      ,C.coeffs_[           5 ].real(),1e-12);
                CHECK_CLOSE(  -80.000000000000000      ,C.coeffs_[           5 ].imag(),1e-12);
                CHECK_CLOSE(   3660.0000000000000      ,C.coeffs_[           6 ].real(),1e-12);
                CHECK_CLOSE(  -160.00000000000000      ,C.coeffs_[           6 ].imag(),1e-12);
                CHECK_CLOSE(   4284.0000000000000      ,C.coeffs_[           7 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[           7 ].imag(),1e-12);
                CHECK_CLOSE(   4908.0000000000000      ,C.coeffs_[           8 ].real(),1e-12);
                CHECK_CLOSE(  -96.000000000000000      ,C.coeffs_[           8 ].imag(),1e-12);
                CHECK_CLOSE(   3812.0000000000000      ,C.coeffs_[           9 ].real(),1e-12);
                CHECK_CLOSE(  -176.00000000000000      ,C.coeffs_[           9 ].imag(),1e-12);
                CHECK_CLOSE(   4468.0000000000000      ,C.coeffs_[          10 ].real(),1e-12);
                CHECK_CLOSE(  -144.00000000000000      ,C.coeffs_[          10 ].imag(),1e-12);
                CHECK_CLOSE(   5124.0000000000000      ,C.coeffs_[          11 ].real(),1e-12);
                CHECK_CLOSE(  -112.00000000000000      ,C.coeffs_[          11 ].imag(),1e-12);
                CHECK_CLOSE(   3964.0000000000000      ,C.coeffs_[          12 ].real(),1e-12);
                CHECK_CLOSE(  -192.00000000000000      ,C.coeffs_[          12 ].imag(),1e-12);
                CHECK_CLOSE(   4652.0000000000000      ,C.coeffs_[          13 ].real(),1e-12);
                CHECK_CLOSE(  -160.00000000000000      ,C.coeffs_[          13 ].imag(),1e-12);
                CHECK_CLOSE(   5340.0000000000000      ,C.coeffs_[          14 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[          14 ].imag(),1e-12);


    }

    TEST(threeResultingIndexContraction){
        const Tensorcd A{1,2,3,4};
        const Tensorcd B{1,5,4,2,3};
        Tensorcd C{1,2,5,4};

        const std::vector<std::pair<int,int>> contraction_pairs{{0,0},{1,3},{3,2}};

        // prepare some demo data
        for(int i = 0; i < A.shape().totalDimension(); ++i){
            A.coeffs_[i] = std::complex<double>(i-1.,i+1.);
        }

        for(int i = 0; i < B.shape().totalDimension(); ++i){
            B.coeffs_[i] = std::complex<double>(i-1.,i+1.);
        }

        general_contraction(A,B,C,contraction_pairs);

        // result should be 3x5x3
                CHECK_EQUAL(3,C.shape().dimensions()[0]);
                CHECK_EQUAL(5,C.shape().dimensions()[1]);
                CHECK_EQUAL(3,C.shape().dimensions()[2]);

                CHECK_CLOSE(   3356.0000000000000      ,C.coeffs_[           0 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[           0 ].imag(),1e-12);
                CHECK_CLOSE(   9436.0000000000000      ,C.coeffs_[          15 ].real(),1e-12);
                CHECK_CLOSE(  -768.00000000000000      ,C.coeffs_[          15 ].imag(),1e-12);
                CHECK_CLOSE(   15516.000000000000      ,C.coeffs_[          30 ].real(),1e-12);
                CHECK_CLOSE(  -1408.0000000000000      ,C.coeffs_[          30 ].imag(),1e-12);
                CHECK_CLOSE(   3916.0000000000000      ,C.coeffs_[           1 ].real(),1e-12);
                CHECK_CLOSE(  -96.000000000000000      ,C.coeffs_[           1 ].imag(),1e-12);
                CHECK_CLOSE(   11276.000000000000      ,C.coeffs_[          16 ].real(),1e-12);
                CHECK_CLOSE(  -736.00000000000000      ,C.coeffs_[          16 ].imag(),1e-12);
                CHECK_CLOSE(   18636.000000000000      ,C.coeffs_[          31 ].real(),1e-12);
                CHECK_CLOSE(  -1376.0000000000000      ,C.coeffs_[          31 ].imag(),1e-12);
                CHECK_CLOSE(   4476.0000000000000      ,C.coeffs_[           2 ].real(),1e-12);
                CHECK_CLOSE(  -64.000000000000000      ,C.coeffs_[           2 ].imag(),1e-12);
                CHECK_CLOSE(   13116.000000000000      ,C.coeffs_[          17 ].real(),1e-12);
                CHECK_CLOSE(  -704.00000000000000      ,C.coeffs_[          17 ].imag(),1e-12);
                CHECK_CLOSE(   21756.000000000000      ,C.coeffs_[          32 ].real(),1e-12);
                CHECK_CLOSE(  -1344.0000000000000      ,C.coeffs_[          32 ].imag(),1e-12);
                CHECK_CLOSE(   3508.0000000000000      ,C.coeffs_[           3 ].real(),1e-12);
                CHECK_CLOSE(  -144.00000000000000      ,C.coeffs_[           3 ].imag(),1e-12);
                CHECK_CLOSE(   9588.0000000000000      ,C.coeffs_[          18 ].real(),1e-12);
                CHECK_CLOSE(  -784.00000000000000      ,C.coeffs_[          18 ].imag(),1e-12);
                CHECK_CLOSE(   15668.000000000000      ,C.coeffs_[          33 ].real(),1e-12);
                CHECK_CLOSE(  -1424.0000000000000      ,C.coeffs_[          33 ].imag(),1e-12);
                CHECK_CLOSE(   4100.0000000000000      ,C.coeffs_[           4 ].real(),1e-12);
                CHECK_CLOSE(  -112.00000000000000      ,C.coeffs_[           4 ].imag(),1e-12);
                CHECK_CLOSE(   11460.000000000000      ,C.coeffs_[          19 ].real(),1e-12);
                CHECK_CLOSE(  -752.00000000000000      ,C.coeffs_[          19 ].imag(),1e-12);
                CHECK_CLOSE(   18820.000000000000      ,C.coeffs_[          34 ].real(),1e-12);
                CHECK_CLOSE(  -1392.0000000000000      ,C.coeffs_[          34 ].imag(),1e-12);
                CHECK_CLOSE(   4692.0000000000000      ,C.coeffs_[           5 ].real(),1e-12);
                CHECK_CLOSE(  -80.000000000000000      ,C.coeffs_[           5 ].imag(),1e-12);
                CHECK_CLOSE(   13332.000000000000      ,C.coeffs_[          20 ].real(),1e-12);
                CHECK_CLOSE(  -720.00000000000000      ,C.coeffs_[          20 ].imag(),1e-12);
                CHECK_CLOSE(   21972.000000000000      ,C.coeffs_[          35 ].real(),1e-12);
                CHECK_CLOSE(  -1360.0000000000000      ,C.coeffs_[          35 ].imag(),1e-12);
                CHECK_CLOSE(   3660.0000000000000      ,C.coeffs_[           6 ].real(),1e-12);
                CHECK_CLOSE(  -160.00000000000000      ,C.coeffs_[           6 ].imag(),1e-12);
                CHECK_CLOSE(   9740.0000000000000      ,C.coeffs_[          21 ].real(),1e-12);
                CHECK_CLOSE(  -800.00000000000000      ,C.coeffs_[          21 ].imag(),1e-12);
                CHECK_CLOSE(   15820.000000000000      ,C.coeffs_[          36 ].real(),1e-12);
                CHECK_CLOSE(  -1440.0000000000000      ,C.coeffs_[          36 ].imag(),1e-12);
                CHECK_CLOSE(   4284.0000000000000      ,C.coeffs_[           7 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[           7 ].imag(),1e-12);
                CHECK_CLOSE(   11644.000000000000      ,C.coeffs_[          22 ].real(),1e-12);
                CHECK_CLOSE(  -768.00000000000000      ,C.coeffs_[          22 ].imag(),1e-12);
                CHECK_CLOSE(   19004.000000000000      ,C.coeffs_[          37 ].real(),1e-12);
                CHECK_CLOSE(  -1408.0000000000000      ,C.coeffs_[          37 ].imag(),1e-12);
                CHECK_CLOSE(   4908.0000000000000      ,C.coeffs_[           8 ].real(),1e-12);
                CHECK_CLOSE(  -96.000000000000000      ,C.coeffs_[           8 ].imag(),1e-12);
                CHECK_CLOSE(   13548.000000000000      ,C.coeffs_[          23 ].real(),1e-12);
                CHECK_CLOSE(  -736.00000000000000      ,C.coeffs_[          23 ].imag(),1e-12);
                CHECK_CLOSE(   22188.000000000000      ,C.coeffs_[          38 ].real(),1e-12);
                CHECK_CLOSE(  -1376.0000000000000      ,C.coeffs_[          38 ].imag(),1e-12);
                CHECK_CLOSE(   3812.0000000000000      ,C.coeffs_[           9 ].real(),1e-12);
                CHECK_CLOSE(  -176.00000000000000      ,C.coeffs_[           9 ].imag(),1e-12);
                CHECK_CLOSE(   9892.0000000000000      ,C.coeffs_[          24 ].real(),1e-12);
                CHECK_CLOSE(  -816.00000000000000      ,C.coeffs_[          24 ].imag(),1e-12);
                CHECK_CLOSE(   15972.000000000000      ,C.coeffs_[          39 ].real(),1e-12);
                CHECK_CLOSE(  -1456.0000000000000      ,C.coeffs_[          39 ].imag(),1e-12);
                CHECK_CLOSE(   4468.0000000000000      ,C.coeffs_[          10 ].real(),1e-12);
                CHECK_CLOSE(  -144.00000000000000      ,C.coeffs_[          10 ].imag(),1e-12);
                CHECK_CLOSE(   11828.000000000000      ,C.coeffs_[          25 ].real(),1e-12);
                CHECK_CLOSE(  -784.00000000000000      ,C.coeffs_[          25 ].imag(),1e-12);
                CHECK_CLOSE(   19188.000000000000      ,C.coeffs_[          40 ].real(),1e-12);
                CHECK_CLOSE(  -1424.0000000000000      ,C.coeffs_[          40 ].imag(),1e-12);
                CHECK_CLOSE(   5124.0000000000000      ,C.coeffs_[          11 ].real(),1e-12);
                CHECK_CLOSE(  -112.00000000000000      ,C.coeffs_[          11 ].imag(),1e-12);
                CHECK_CLOSE(   13764.000000000000      ,C.coeffs_[          26 ].real(),1e-12);
                CHECK_CLOSE(  -752.00000000000000      ,C.coeffs_[          26 ].imag(),1e-12);
                CHECK_CLOSE(   22404.000000000000      ,C.coeffs_[          41 ].real(),1e-12);
                CHECK_CLOSE(  -1392.0000000000000      ,C.coeffs_[          41 ].imag(),1e-12);
                CHECK_CLOSE(   3964.0000000000000      ,C.coeffs_[          12 ].real(),1e-12);
                CHECK_CLOSE(  -192.00000000000000      ,C.coeffs_[          12 ].imag(),1e-12);
                CHECK_CLOSE(   10044.000000000000      ,C.coeffs_[          27 ].real(),1e-12);
                CHECK_CLOSE(  -832.00000000000000      ,C.coeffs_[          27 ].imag(),1e-12);
                CHECK_CLOSE(   16124.000000000000      ,C.coeffs_[          42 ].real(),1e-12);
                CHECK_CLOSE(  -1472.0000000000000      ,C.coeffs_[          42 ].imag(),1e-12);
                CHECK_CLOSE(   4652.0000000000000      ,C.coeffs_[          13 ].real(),1e-12);
                CHECK_CLOSE(  -160.00000000000000      ,C.coeffs_[          13 ].imag(),1e-12);
                CHECK_CLOSE(   12012.000000000000      ,C.coeffs_[          28 ].real(),1e-12);
                CHECK_CLOSE(  -800.00000000000000      ,C.coeffs_[          28 ].imag(),1e-12);
                CHECK_CLOSE(   19372.000000000000      ,C.coeffs_[          43 ].real(),1e-12);
                CHECK_CLOSE(  -1440.0000000000000      ,C.coeffs_[          43 ].imag(),1e-12);
                CHECK_CLOSE(   5340.0000000000000      ,C.coeffs_[          14 ].real(),1e-12);
                CHECK_CLOSE(  -128.00000000000000      ,C.coeffs_[          14 ].imag(),1e-12);
                CHECK_CLOSE(   13980.000000000000      ,C.coeffs_[          29 ].real(),1e-12);
                CHECK_CLOSE(  -768.00000000000000      ,C.coeffs_[          29 ].imag(),1e-12);
                CHECK_CLOSE(   22620.000000000000      ,C.coeffs_[          44 ].real(),1e-12);
                CHECK_CLOSE(  -1408.0000000000000      ,C.coeffs_[          44 ].imag(),1e-12);



    }

}

