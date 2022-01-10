//
// Created by Roman Ellerbrock on 12/15/21.
//
#include "Benchmark.h"

void benchmark(const std::function<void()>& f, const std::string& fName, size_t nSample) {

	Timings t;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	for (size_t i = 0; i < nSample; ++i) {
		start = std::chrono::system_clock::now();
		f();
		end = std::chrono::system_clock::now();
		t.emplace_back(std::chrono::duration_cast<std::chrono::microseconds>(end - start).count());
	}
	t.report(fName);
}

