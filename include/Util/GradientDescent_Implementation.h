//
// Created by Roman Ellerbrock on 2/24/20.
//

#include "Util/GradientDescent.h"

double lr_schedule(double learning_rate, double time, double max_time) {
	return learning_rate;
}

template<class Model, class Parameters>
void GradientDescent(Model& model, double learning_rate, size_t num_iter) {
	for (size_t iter = 0; iter < num_iter; ++iter) {
		learning_rate = lr_schedule(learning_rate, iter, num_iter);
		model.print();
		Parameters grad = model.gradient();
		model.parameters() -= grad * learning_rate;
		cout << "Gradient:";
		grad.print();
	}
}



