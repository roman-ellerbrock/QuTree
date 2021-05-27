//
// Created by Roman Ellerbrock on 5/22/21.
//
#include "benchmark_tensor.h"
#include "benchmark_helper.h"
#include "Core/TensorBLAS.h"

namespace benchmark {
	auto sampleTranspose(Matrixcd& dest, const Matrixcd& src,
		size_t nsample, size_t dim1, size_t dim2, size_t blocksize) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			transpose2<complex<double>,4>(&dest[0], &src[0], dim1, dim2);
//			transpose(&dest[0], &src[0], dim1, dim2);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto sampleTransposeAB(Tensorcd& dest, const Tensorcd& src,
		size_t nsample, size_t A, size_t B, size_t C, size_t blocksize) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			transposeAB(&dest[0], &src[0], A, B, C);
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto sampleMatrixTensor(Tensorcd& B, const Matrixcd& S, const Tensorcd& A,
		size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
		Tensorcd B2 = B;
		Tensorcd A2(A);
		for (size_t n = 0; n < nsample; ++n) {
			B.zero();
//			B2.zero();
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
//			matrixTensor1(B, S, A, bef, act, act, aft, true);
			matrixTensor2(B, S, A, B2, bef, act, act, aft, true);
//			matrixTensor3(B, S, A, A2, B2, bef, act, act, aft, true);
//			matrixTensor(B2, S, A, bef, act, act, aft, true);
//			cout << "res = " << residual(B, B2) << endl;
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto sampleTensorHoleProduct(Matrixcd& S, const Tensorcd& B, const Tensorcd& A,
		size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
		Tensorcd workA = A;
		Tensorcd workB = B;
		for (size_t n = 0; n < nsample; ++n) {
			S.zero();
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			contraction2(S, A, B, workA, workB, bef, act, act, aft, true);
//			auto S2 = contraction(A, B, 0);
//			contraction(S, A, B, 2);
//			cout << residual(S, S2) << endl;
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	void screenTranspose(mt19937& gen, ostream& os, size_t nsample) {
		size_t blocksize = 8;
		size_t chunk = 32;
		size_t max = 512;
		for (size_t dim = chunk; dim <= max; dim += chunk) {
			Matrixcd src(dim, dim);
			Matrixcd dest(dim, dim);
			Tensor_Extension::generate(src, gen);
			auto stat = sampleTranspose(dest, src, nsample, dim, dim, blocksize);
			os << dim << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

	void screenDimensionTransposeAB(mt19937& gen, ostream& os, size_t nsample) {
		size_t order = 3;
		size_t mode = 1;
		size_t blocksize = 4;
		for (size_t dim = 50; dim <= 250; dim += 20) {
			auto tdim = make_TensorDim(order, dim);
			Tensorcd A(tdim, false);
			Tensorcd B(tdim, false);
			Matrixcd S(dim, dim);
			Tensor_Extension::generate(A, gen);
			Tensor_Extension::generate(B, gen);
			size_t aft = tdim.after(mode);
			size_t act = tdim[mode];
			size_t bef = tdim.before(mode);

			auto stat = sampleTransposeAB(A, B, nsample, bef, act, aft, blocksize);
			os << dim << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

	void screenDimensionMatrixTensor(mt19937& gen, ostream& os, size_t nsample) {
		size_t order = 5;
//		size_t dim = 10;
		size_t mode = 4;
//		for (size_t dim = 50; dim <= 250; dim += 20) { // order 3
		for (size_t dim = 8; dim <= 36; dim += 4) { // order 5
			auto tdim = make_TensorDim(order, dim);
			Tensorcd A(tdim, false);
			Matrixcd S(dim, dim);
			Tensor_Extension::generate(A, gen);
			Tensor_Extension::generate(S, gen);
			Tensorcd B(tdim, true);
			size_t aft = tdim.after(mode);
			size_t act = tdim[mode];
			size_t bef = tdim.before(mode);

			auto stat = sampleMatrixTensor(B, S, A, nsample, bef, act, aft);
			os << dim << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

	void screenDimensionTensorHoleProduct(mt19937& gen, ostream& os, size_t nsample) {
		size_t order = 5;
		size_t mode = 2;
//		for (size_t dim = 50; dim <= 250; dim += 20) {
		for (size_t dim = 8; dim <= 48; dim += 4) { // order 5
			auto tdim = make_TensorDim(order, dim);
			Tensorcd A(tdim, false);
			Tensorcd B(tdim, false);
			Matrixcd S(dim, dim);
			Tensor_Extension::generate(A, gen);
			Tensor_Extension::generate(B, gen);
			size_t aft = tdim.after(mode);
			size_t act = tdim[mode];
			size_t bef = tdim.before(mode);

			auto stat = sampleTensorHoleProduct(S, A, B, nsample, bef, act, aft);
			os << dim << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

}
