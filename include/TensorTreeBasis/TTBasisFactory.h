//
// Created by Roman Ellerbrock on 2/2/20.
//

#ifndef TTBASISFACTORY_H
#define TTBASISFACTORY_H
#include "TensorTreeBasis/TensorTreeBasis.h"

namespace TTBasisFactory {
	TTBasis TensorTrain(size_t nLeaves, size_t dimLeaves, size_t dimNodes, size_t leafType);
}

#endif //TTBASISFACTORY_H