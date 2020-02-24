//
// Created by Roman Ellerbrock on 2/24/20.
//

#include "UnitTest++/UnitTest++.h"
#include "Util/GradientDescent_Implementation.h"
#include "Core/Vector.h"

SUITE(GradientDescent) {
	class Interface : public Vectord {
		public:

		Interface(): Vectord(2) {
			operator()(0) = 0.5;
			operator()(1) = -0.2;
			alpha_ = 0.5;
			beta_ = 0.4;
		}

		~Interface() = default;

		double evaluate(const Vectord& x) const {
			return -exp(-alpha_ * x(0) * x(0) - beta_ * x(1) * x(1)) + 1.;
		}

		Vectord Gradient() const {
			Vectord grad(2);
			grad(0) = 2 * alpha_ * operator()(0) * evaluate(*this);
			grad(1) = 2 * beta_ * operator()(1) * evaluate(*this);
			return grad;
		}

		void print()const {
			cout << "x = ";
			Vectord::print();
			cout << "f(x) = " << evaluate(*this) << endl;
		}

		private:
		double alpha_;
		double beta_;
	};

	TEST(GradientDescent1) {
		Interface myInterface;
		myInterface.print();
		GradientDescent<Interface, Vectord>(myInterface, 1., 100);
		myInterface.print();

	}
}