//
// Created by Roman Ellerbrock on 2/3/20.
//
#include "../tests/benchmarks/benchmark_Tensor.h"

namespace benchmark_tensor {
	TensorDim make_TensorDim(size_t order, size_t dim) {
		assert(order > 1);
		assert(dim > 0);
		vector<size_t> dims;
		for (size_t i = 0; i < order - 1; ++i) {
			dims.emplace_back(dim);
		}
		return TensorDim(dims, dim);
	}

	void dot_product(mt19937& gen, size_t dim, size_t max_order, ostream& os) {
		assert(max_order >= 2);
		for (size_t order = 2; order <= max_order; ++order) {
			auto tdim =  make_TensorDim(order, dim);
			/// Initialize memory
			std::chrono::time_point<std::chrono::system_clock> allocate_start, allocate_end;
			allocate_start = std::chrono::system_clock::now();
			Tensorcd A(tdim, false);
			Tensorcd B(tdim, false);
			Tensor_Extension::Generate(A, gen);
			Tensor_Extension::Generate(B, gen);
			Matrixcd S(dim, dim);
			size_t bef = tdim.GetDimPart();
			size_t act = tdim.GetNumTensor();

			allocate_end = std::chrono::system_clock::now();
			std::chrono::duration<double> allocate_elapsed_seconds = allocate_end - allocate_start;
			os << "Preparation time: " << allocate_elapsed_seconds.count() << endl;
			std::chrono::time_point<std::chrono::system_clock> start, end;
			/// Run hole-product
			start = std::chrono::system_clock::now();
			TensorHoleProduct(S, A, B, bef, act, act, 1);
			end = std::chrono::system_clock::now();

			std::chrono::duration<double> elapsed_seconds = end - start;
			os << dim << "\t" << order << "\t" << elapsed_seconds.count() << endl;
		}
	}

	void hole_product(mt19937& gen, size_t dim, size_t max_order, ostream& os) {
		assert(max_order >= 2);
		for (size_t order = 3; order <= max_order; order += 2) {
			os << dim << "\t" << order;
			size_t step = (order)/2;

			step = (step == 0) ? 1 : step;
			for (size_t mode = 0; mode < order; mode += step) {
				auto tdim = make_TensorDim(order, dim);
				/// Initialize memory
				std::chrono::time_point<std::chrono::system_clock> allocate_start, allocate_end;
				allocate_start = std::chrono::system_clock::now();
				Tensorcd A(tdim, false);
				Tensorcd B(tdim, false);
				Tensor_Extension::Generate(A, gen);
				Tensor_Extension::Generate(B, gen);
				Matrixcd S(dim, dim);
				size_t aft = tdim.After(mode);
				size_t act = tdim.Active(mode);
				size_t bef = tdim.Before(mode);

				allocate_end = std::chrono::system_clock::now();
				std::chrono::duration<double> allocate_elapsed_seconds = allocate_end - allocate_start;
//				os << "Preparation time: " << allocate_elapsed_seconds.count() << endl;

				/// Run hole-product
				std::chrono::time_point<std::chrono::system_clock> start, end;
				start = std::chrono::system_clock::now();
				TensorHoleProduct(S, A, B, bef, act, act, aft);
				end = std::chrono::system_clock::now();

				auto duration = chrono::duration_cast<chrono::microseconds>(end-start).count();
				os << "\t" <<  duration;
			}
			os << endl;
		}
	}

	void matrix_tensor(mt19937& gen, size_t dim, size_t max_order, ostream& os) {
		assert(max_order >= 2);
		for (size_t order = 3; order <= max_order; order += 2) {
			size_t step = (order)/2;

			os << dim << "\t" << order;
			step = (step == 0) ? 1 : step;
			for (size_t mode = 0; mode < order; mode += step) {
				/// Initialize memory
				std::chrono::time_point<std::chrono::system_clock> allocate_start, allocate_end;
				allocate_start = std::chrono::system_clock::now();
				auto tdim = make_TensorDim(order, dim);
				Tensorcd A(tdim, false);
				Matrixcd S(dim, dim);
				Tensor_Extension::Generate(A, gen);
				Tensor_Extension::Generate(S, gen);
				Tensorcd B(tdim, true);

				auto start = std::chrono::system_clock::now();
				multAB(B, S, A, mode);
				auto end = std::chrono::system_clock::now();
				auto duration = chrono::duration_cast<chrono::microseconds>(end-start).count();
				os << "\t" <<  duration;
			}
			os << endl;
		}
	}

	void run() {
		mt19937 gen(1989);
/*		dot_product(gen, 2, 25);
		dot_product(gen, 32, 5);
		dot_product(gen, 64, 4);
*/

		size_t dim = 2;
		size_t max_order =  29;
		hole_product(gen, dim, max_order, cout);

		matrix_tensor(gen, dim, max_order, cout);

	}
}

