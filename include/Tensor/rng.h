//
// Created by Roman Ellerbrock on 11/19/21.
//

#ifndef RNG_H
#define RNG_H

#include <random>

namespace rng {
	static uniform_real_distribution uniform(-1., 1.);
	static uniform_real_distribution normal(-1., 1.);
	static std::mt19937 gen(time(nullptr));
	static std::mt19937 genConstSeed(523838);
}

#endif //RNG_H
