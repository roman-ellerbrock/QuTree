//
// Created by Roman Ellerbrock on 2/4/20.
//
#include "benchmark_helper.h"
#include <blas/flops.hh>
#include "Tensor/TensorBLAS2.h"

namespace benchmark {
	TensorShape make_TensorDim(size_t order, size_t dim) {
		assert(order > 1);
		assert(dim > 0);
		vector<size_t> dims;
		for (size_t i = 0; i < order; ++i) {
			dims.emplace_back(dim);
		}
		return TensorShape(dims);
	}

	template <typename T>
	double calc_gflops(const TensorShape& dim, size_t k) {
		double inactive = dim.before(k) * dim.after(k);
		double n = dim[k] * dim[k] * inactive;
		auto mul_ops = blas::FlopTraits<T>::mul_ops;
		auto add_ops = blas::FlopTraits<T>::add_ops;
		return 1e-9 * (mul_ops * n + add_ops * n);
	}
	template double calc_gflops<double>(const TensorShape& dim, size_t k);
	template double calc_gflops<complex<double>>(const TensorShape& dim, size_t k);

	template <typename T>
	double calc_mem(const TensorShape& dim, size_t k) {
		double inactive = dim.before(k) * dim.after(k);
		double n = dim[k] * dim[k] * inactive;
		auto mul_ops = blas::FlopTraits<T>::mul_ops;
		auto add_ops = blas::FlopTraits<T>::add_ops;
		/// read/write h, read A & B
		return 1e-9 * (2 * dim[k] * dim[k] + 2 * dim.totalDimension()) * sizeof(T);
	}
	template double calc_mem<double>(const TensorShape& dim, size_t k);
	template double calc_mem<complex<double>>(const TensorShape& dim, size_t k);


	template <typename T>
	double calc_mem_matrixTensor(const TensorShape& dim, size_t k) {
		double inactive = dim.before(k) * dim.after(k);
		double n = dim[k] * dim[k] * inactive;
		auto mul_ops = blas::FlopTraits<T>::mul_ops;
		auto add_ops = blas::FlopTraits<T>::add_ops;
		/// B = h * A
		/// read/write B, read h & A
		return 1e-9 * (dim[k] * dim[k] + 3 * dim.totalDimension()) * sizeof(T);
	}

	template double calc_mem_matrixTensor<double>(const TensorShape& dim, size_t k);
	template double calc_mem_matrixTensor<complex<double>>(const TensorShape& dim, size_t k);


	pair<double, double> statistic_helper(const vector<chrono::microseconds>& dur) {
		double mean = 0;
		for (auto x : dur) {
			mean += x.count();
		}
		mean /= (double) dur.size();

		double std_dev = 0;
		for (auto x : dur) {
			std_dev += pow(x.count() - mean, 2);
		}
		std_dev = sqrt(std_dev / ((double) dur.size() - 1));
		return pair<double, double>(mean, std_dev);
	}

	pair<double, double> statistic_helper(const vector<double>& dur) {
		double mean = 0;
		for (auto x : dur) {
			mean += x;
		}
		mean /= (double) dur.size();

		double std_dev = 0;
		for (auto x : dur) {
			std_dev += pow(x - mean, 2);
		}
		std_dev = sqrt(std_dev / ((double) dur.size() - 1));
		return pair<double, double>(mean, std_dev);
	}

}