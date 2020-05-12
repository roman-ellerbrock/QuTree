//
// Created by Roman Ellerbrock on 3/3/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Util/RungeKutta4.h"
#include "Core/Vector.h"
#include "Util/QMConstants.h"

SUITE(Integrators) {
	class HOInterface {
	public:
		HOInterface() = default;
		~HOInterface() = default;

		Vectord Derivative(double t, const Vectord& y) {
			Vectord dy(2);
			dy(0) = -y(1); // dp
			dy(1) = y(0); // dq
			return dy;
		}

		void summary(double t, const Vectord& y) {
/*			cout << t << "\t";
			y.print();*/
		}
	};

	TEST(RK4) {
		/// Integrate a 1d Harmonic Oscillator for one circle
		Vectord y(2);
		y(0) = 0.;
		y(1) = -1.;
		double t = 0.;
		double t_end = QM::two_pi;
		double h = 0.0005;
		HOInterface I;

		auto y_correct = y;
		RungeKutta4::Integrate<HOInterface, Vectord, double>(t, t_end, h, y, I);
		auto residual = Residual(y_correct, y);
		CHECK_CLOSE(0., residual, 10. * h);
	}

}