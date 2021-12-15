//
// Created by Roman Ellerbrock on 12/15/21.
//

#ifndef BENCHMARKER_H
#define BENCHMARKER_H


class Benchmarker {
public:
	Benchmarker() = default;
	~Benchmarker() = default;


	double mean() const {
		double mean = 0;
		for (auto x : dur_) {
			mean += x.count();
		}
		mean /= (double) dur_.size();

		return mean;
	}

	double stdDev() const {
		double std_dev = 0;
		double m = mean();
		for (auto x : dur) {
			std_dev += pow(x.count() - m, 2);
		}
		std_dev = sqrt(std_dev / ((double) dur.size() - 1));
		return std_dev;
	}
	vector<chrono::microseconds>& dur_;
};


#endif //BENCHMARKER_H
