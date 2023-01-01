//
// Created by Roman Ellerbrock on 12/31/22.
//

#ifndef QUTREE_RNG_H
#define QUTREE_RNG_H
#include <random>

namespace qutree {

	static std::mt19937 rng;

	void init_rng() {
		rng = std::mt19937(time(NULL));
	}
}

#endif //QUTREE_RNG_H
