//
// Created by Roman Ellerbrock on 3/3/20.
//

#ifndef RUNGEKUTTA4_H
#define RUNGEKUTTA4_H
#include <iostream>
#include <fstream>
#include <functional>

namespace RungeKutta4 {
	using namespace std;

	/**
	 * \brief This class performs a 4th order Runge-Kutta integration for a given Model, Container and Type.
	 * @tparam Container is a class holding all required memory to obtain the derivatives
	 * @tparam Model is the object you want to integrate. E.g. a Tensor, a TensorTree, etc.
	 * @tparam U is the underlying datatype, E.g. double, complex<double>, float, etc.
	 */
	template<class Container, class Model, typename U>
	class RK_integrator {
	public:
		RK_integrator() = default;
		RK_integrator(size_t dim, const Model& empty) : k1(empty), k2(empty), k3(empty), k4(empty) {
		}
		~RK_integrator() = default;

		void Integrate(Model& y, double& t, double tend, double& dt, double eps,
			function<void(Container&, double, Model&, const Model&)> ddx,
			function<double(Container&, Model&, Model&)> err, Container& I) {

			while (t + 1e-7 < tend) {
				dt = std::min(dt, tend - t);
				step(t, y, dt, I, ddx);
			}
		}

/*		void integrate(double t, double t_end, double dt,
			Model& y, Container& container) {
			while (t + 1e-7 < t_end) {
				dt = std::min(dt, t_end - t);
				step(t, y, dt, container);
			}
		}*/

		/**
		 *
		 * @param t is current time
		 * @param y is current Model (Tensor, TensorTree, ...)
		 * @param dt is stepsize in time
		 * @param I is the memory container
		 * @param ddx
		 * @param container
		 */
		void step(U& t, Model& y, U dt, Container& I,
			function<void(Container&, double, Model&, const Model&)> ddx) {
			// Runge kutta values
//			Model k1 = I.Derivative(t, y) * h;
			ddx(I, t, k1, y);
			k1 *= dt;

//			Model k2 = I.Derivative(t + h / 2., y + k1 / 2.) * h;
			ddx(I, t + dt / 2., k2, y + (0.5 * k1));
			k2 *= dt;

//			Model k3 = I.Derivative(t + h / 2., y + k2 / 2.) * h;
			ddx(I, t + dt / 2., k3, y + (0.5 * k2));
			k3 *= dt;

//			Model k4 = I.Derivative(t + h, y + k3) * h;
			ddx(I, t + dt, k4, y + k3);
			k4 *= dt;

			/// add all together to obtain the result
			// y += 1./6. * (k1 + 2. * k2 + 2. * k3 + k4);
			y +=  1. / 6. * (k1 + 2. * k2 + 2. * k3 + k4);
			t += dt;
		}


		Model k1, k2, k3, k4;
	};

	/*
	template<class Interface, class Model, typename T>
	void step(T& t, Model& y, T h, Interface& I) {
		// Runge kutta values
		Model k1 = I.Derivative(t, y) * h;
		Model k2 = I.Derivative(t + h / 2., y + k1 / 2.) * h;
		Model k3 = I.Derivative(t + h / 2., y + k2 / 2.) * h;
		Model k4 = I.Derivative(t + h, y + k3) * h;

		y += 1./6. * (k1 + 2. * k2 + 2. * k3 + k4);
		t += h;
	}

	template<class Interface, class Model, typename T>
	void integrate(double t, double t_end, double h,
		Model& y, Interface& I) {
//		std::ofstream os("ho.dat");
		while (t + 1e-7 < t_end) {
			h = std::min(h, t_end - t);
			step(t, y, h, I);
//			I.summary(t, y);
		}
	}
*/
}

#endif //RUNGEKUTTA4_H
