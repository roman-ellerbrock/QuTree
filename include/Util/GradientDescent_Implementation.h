//
// Created by Roman Ellerbrock on 2/24/20.
//

#include "Util/GradientDescent.h"

double lr_schedule(double learning_rate, double time, double max_time) {
	return learning_rate;
}

template<class LearningInterface, class Func>
void GradientDescent(LearningInterface& func, double learning_rate, size_t num_iter) {
	for (size_t iter = 0; iter < num_iter; ++iter) {
		learning_rate = lr_schedule(learning_rate, iter, num_iter);
		func.print();
		Func grad = func.Gradient();
		func -= grad * learning_rate;
		cout << "Gradient:";
		grad.print();
	}
}



