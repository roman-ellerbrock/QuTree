//
// Created by Roman Ellerbrock on 2/4/20.
//
#include "benchmark_helper.h"

namespace benchmark {
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