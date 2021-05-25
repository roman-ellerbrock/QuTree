//
// Created by Roman Ellerbrock on 5/22/21.
//
#include "benchmark_tensor.h"
#include "benchmark_helper.h"
#include "Core/TensorMath.h"

namespace benchmark {
	auto matrix_tensor1(Tensorcd& B, const Matrixcd& S, const Tensorcd& A,
		size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
//		Tensorcd B2 = B;
		for (size_t n = 0; n < nsample; ++n) {
			B.zero();
//			B2.zero();
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			matrixTensor(B, S, A, bef, act, act, aft, false);
//			matrixTensor(B2, S, A, bef, act, act, aft, false);
//			cout << "res = " << residual(B, B2) << endl;
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	auto tensor_tensor1(Matrixcd& S, const Tensorcd& B,  const Tensorcd& A,
		size_t nsample, size_t bef, size_t act, size_t aft) {
		vector<chrono::microseconds> duration_vec;
		for (size_t n = 0; n < nsample; ++n) {
			S.zero();
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			contraction2(S, A, B, bef, act, act, aft, false);
//			auto S2 = contraction(A, B, 0);
//			cout << residual(S, S2) << endl;
			end = std::chrono::system_clock::now();
			duration_vec.emplace_back(chrono::duration_cast<chrono::microseconds>(end - start).count());
		}

		return statistic_helper(duration_vec);
	}

	void screen_matrixtensor_optimization(mt19937& gen, ostream& os, size_t nsample) {
		size_t order = 3;
		size_t dim = 10;
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A(tdim, false);
		Matrixcd S(dim, dim);
		Tensor_Extension::generate(A, gen);
		Tensor_Extension::generate(S, gen);
		Tensorcd B(tdim, true);
		for (size_t mode = 0; mode < order; ++mode) {
			size_t aft = tdim.after(mode);
			size_t act = tdim[mode];
			size_t bef = tdim.before(mode);

			auto stat = matrix_tensor1(B, S, A, nsample, bef, act, aft);
			os << "mode = " << mode << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

	void screen_contraction_optimization(mt19937& gen, ostream& os, size_t nsample) {
		size_t order = 3;
		size_t dim = 10;
		auto tdim = make_TensorDim(order, dim);
		Tensorcd A(tdim, false);
		Tensorcd B(tdim, false);
		Matrixcd S(dim, dim);
		Tensor_Extension::generate(A, gen);
		Tensor_Extension::generate(B, gen);
		for (size_t mode = 0; mode < order; ++mode) {
			size_t aft = tdim.after(mode);
			size_t act = tdim[mode];
			size_t bef = tdim.before(mode);

			auto stat = tensor_tensor1(S, A, B, nsample, bef, act, aft);
			os << "mode = " << mode << "\t" << stat.first / 1000. << "\t" << stat.second / 1000. << endl;
		}
	}

}
