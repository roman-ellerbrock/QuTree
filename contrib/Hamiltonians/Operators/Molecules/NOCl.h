//
// Created by Roman Ellerbrock on 7/1/20.
//
#pragma once
#include "Core/stdafx.h"
#include "Core/Vector.h"
#include "TreeOperators/LeafFunction.h"
#include "TreeOperators/SumOfProductsOperator.h"

namespace Operator {
	SOPcd NOCl_KEO();
	SOPcd NOCl_V();
	SOPcd NOCl_H(bool potential);
}
