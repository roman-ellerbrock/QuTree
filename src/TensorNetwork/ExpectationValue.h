
#ifndef EXPECTATIONVALUE_H
#define EXPECTATIONVALUE_H
#include "TensorNetwork/contractions.h"

Tensorcd expectationValue(const TensorTreecd& Psi, const SOPcd& H, const Tree& tree);

#endif // EXPECTATIONVALUE_H