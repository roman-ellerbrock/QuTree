//
// Created by Roman Ellerbrock on 12/7/22.
//

#ifndef ACTIVECOUNTER_H
#define ACTIVECOUNTER_H
#include "TreeClasses/NodeAttribute.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeClasses/SparseMatrixTree.h"

size_t nActives(const SparseMatrixTreecd& hmat, const MLOcd& M, const SparseMatrixTreecd& hcon,
	const Node& node);

/**
 * \brief Counts how many nodes are active per summand in a SOP
 */
class ActiveCounter : public NodeAttribute<vector<size_t>> {
public:
	ActiveCounter() = default;
	~ActiveCounter() = default;

	void initialize(SparseMatrixTreescd& hmats, SparseMatrixTreescd& hcons,
		const SOPcd& S, const Tree& tree) {

		attributes_.clear();
		for (const Node& node : tree) {
			attributes_.emplace_back(vector<size_t>(S.size()));

			for (size_t l = 0; l < S.size(); ++l) {
				(*this)[node][l] = nActives(hmats[l], S[l], hcons[l], node);
			}
		}

	}
};


#endif //ACTIVECOUNTER_H
