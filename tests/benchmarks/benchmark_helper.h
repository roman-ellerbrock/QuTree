//
// Created by Roman Ellerbrock on 2/4/20.
//

#ifndef BENCHMARK_HELPER_H
#define BENCHMARK_HELPER_H
#include <chrono>
#include "TreeHandling/Tree.h"

namespace benchmark {
	pair<double, double> statistic_helper(const vector<chrono::microseconds>& dur);
}

#endif //BENCHMARK_HELPER_H
