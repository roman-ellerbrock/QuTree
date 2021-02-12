//
// Created by Roman Ellerbrock on 2/24/20.
//

#ifndef GRADIENTDESCENT_H
#define GRADIENTDESCENT_H
#include "Core/stdafx.h"

template<class Interface, class Func>
void gradientDescent(Interface& func, double learning_rate, size_t num_iter);

#endif //GRADIENTDESCENT_H
