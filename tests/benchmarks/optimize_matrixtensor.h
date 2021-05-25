//
// Created by Roman Ellerbrock on 5/22/21.
//

#ifndef OPTIMIZE_MATRIXTENSOR_H
#define OPTIMIZE_MATRIXTENSOR_H

namespace benchmark {
	void screen_matrixtensor_optimization(mt19937& gen, ostream& os, size_t nsample);
	void screen_contraction_optimization(mt19937& gen, ostream& os, size_t nsample);
}

#endif //OPTIMIZE_MATRIXTENSOR_H
