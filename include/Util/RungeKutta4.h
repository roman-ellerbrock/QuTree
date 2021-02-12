//
// Created by Roman Ellerbrock on 3/3/20.
//

#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H
#include <iostream>
#include <fstream>

namespace RungeKutta4 {

	template<class Interface, class Model, typename T>
	void step(T& t, Model& y, T h, Interface& I) {
		// Runge kutta values
		Model k1 = I.Derivative(t, y) * h;
		Model k2 = I.Derivative(t + h / 2., y + k1 / 2.) * h;
		Model k3 = I.Derivative(t + h / 2., y + k2 / 2.) * h;
		Model k4 = I.Derivative(t + h, y + k3) * h;

		y += (k1 + k2 * 2. + k3 * 2. + k4) / 6.;
		t += h;
	}

	template<class Interface, class Model, typename T>
	void integrate(double t, double t_end, double h,
		Model& y, Interface& I) {
		std::ofstream os("ho.dat");
		while (t + 1e-7 < t_end) {
			h = std::min(h, t_end - t);
			step(t, y, h, I);
			I.summary(t, y);
		}
	}

}

#endif //RUNGEKUTTA4_H
