//
// Created by Roman Ellerbrock on 12/15/21.
//

#ifndef BENCHMARK_H
#define BENCHMARK_H
#include <chrono>
#include <vector>
#include <math.h>
#include <string>
#include <iostream>

class Timings : public std::vector<std::chrono::microseconds> {
public:

	[[nodiscard]] double mean() const {
		double m = 0.;
		for (auto x : *this) {
			m += x.count();
		}
		m /= (double) size();
		return m;
	}

	[[nodiscard]] double stdDev() const {
		if (empty()) {
			return 0.;
		}
		double std_dev = 0;
		double m = mean();
		for (auto x : *this) {
			std_dev += pow(x.count() - m, 2);
		}

		std_dev = sqrt(std_dev / ((double) size() - 1));
		return std_dev;
	}

	void report(const std::string& fName) const {
		std::cout << "function: " << fName << " @" << size();
		std::cout << ": <t>=" << mean()/1000. << "+/-";
		std::cout << stdDev()/1000. << " ms" << std::endl;
	}

};

void benchmark(const std::function<void()>& f, const std::string& fName = "", size_t nSample = 10);

#endif //BENCHMARK_H
