//
// Created by Roman Ellerbrock on 2019-07-24.
//

#ifndef QMCONSTANTS_H
#define QMCONSTANTS_H
#include "stdafx.h"

namespace QM {
	// Math constants
	constexpr double pi = 3.1415926535897932384626433832795;
	constexpr double two_pi = 2. * pi;

	constexpr complex<double> im (0., 1.);

	// Energy
	constexpr double cm = 219474.6313705;
	constexpr double eV = 27.2114;
	constexpr double fs = 41.362;

	// imaginary time to Kelvin
	constexpr double Ha_K = 315777.09;
}

#endif //QMCONSTANTS_H
