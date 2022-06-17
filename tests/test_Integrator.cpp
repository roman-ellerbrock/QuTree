//
// Created by Roman Ellerbrock on 3/3/20.
//

#include <gtest/gtest.h>
#include "Util/RungeKutta4.h"
#include "Core/Vector.h"
#include "Util/QMConstants.h"

class HOInterface {
public:
    HOInterface() = default;
    ~HOInterface() = default;

/*		Vectord Derivative(double t, const Vectord& y) {
			Vectord dy(2);
			dy(0) = -y(1); // dp
			dy(1) = y(0); // dq
			return dy;
		}*/

    void derivative(double t, Vectord& dy, const Vectord& y) {
        dy(0) = -y(1); // dp
        dy(1) = y(0); // dq
    }

    double error(Vectord& x, Vectord& y) {
        auto diff = x - y;
        return diff.norm();
    }

    void summary(double t, const Vectord& y) {
/*			cout << t << "\t";
			y.print();*/
    }
};


TEST (Integrators, RK4) {
    /// integrate a 1d Harmonic Oscillator for one circle
    Vectord y(2);
    y(0) = 0.;
    y(1) = -1.;
    double t = 0.;
    double t_end = QM::two_pi;
    double h = 0.005;
    HOInterface I;
    auto y_correct = y;

    RungeKutta4::RK_integrator<HOInterface, Vectord, double> rk(y.dim(), y);
    function<void(HOInterface&, double, Vectord&, const Vectord&)> ddx = &HOInterface::derivative;;
    function<double(HOInterface&, Vectord&, Vectord&)> err = &HOInterface::error;
    double eps = 1e-5; // unused
    rk.Integrate(y, t, t_end, h, eps, ddx, err, I);

    auto res = residual(y_correct, y);
    ASSERT_NEAR(0., res, 10. * h);
}
