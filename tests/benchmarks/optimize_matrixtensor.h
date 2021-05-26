//
// Created by Roman Ellerbrock on 5/22/21.
//

#ifndef OPTIMIZE_MATRIXTENSOR_H
#define OPTIMIZE_MATRIXTENSOR_H

namespace benchmark {
	void screenDimensionMatrixTensor(mt19937& gen, ostream& os, size_t nsample);
	void screenDimensionTensorHoleProduct(mt19937& gen, ostream& os, size_t nsample);
	void screenTranspose(mt19937& gen, ostream& os, size_t nsample);
	void screenDimensionTransposeAB(mt19937& gen, ostream& os, size_t nsample);
}

#endif //OPTIMIZE_MATRIXTENSOR_H
